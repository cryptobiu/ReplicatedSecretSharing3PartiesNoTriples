#ifndef PROTOCOLPARTY_H_
#define PROTOCOLPARTY_H_

#include <stdlib.h>
#include <libscapi/include/primitives/Matrix.hpp>
#include <libscapi/include/circuits/ArithmeticCircuit.hpp>
#include <libscapi/include/comm/MPCCommunication.hpp>
#include <libscapi/include/cryptoInfra/Protocol.hpp>

#include <libscapi/include/primitives/Mersenne.hpp>
//#include "ProtocolTimer.h"
#include <libscapi/include/infra/Measurement.hpp>
#include "HashEncrypt.h"
#include <libscapi/include/infra/Common.hpp>
#include <thread>

#define flag_print false
#define flag_print_timings true
#define flag_print_output true

#define N 3

using namespace std;
using namespace std::chrono;

template <class FieldType>
class ProtocolParty : public Protocol, public HonestMajority, public ThreeParty{
private:
    /**
     * N - number of parties
     * T - number of malicious
     * M - number of gates
     */
    Measurement* timer;
    int currentCirciutLayer = 0;
    int M, m_partyId;
    int numOfInputGates, numOfOutputGates;
    int times; //number of times to run the run function
    int iteration; //number of the current iteration
    string inputsFile, outputFile;
    FieldType inv_3;
    int LEFT, RIGHT;
    int ZERO_PARTY, ONE_PARTY,TWO_PARTY;
    FieldType num_2; // number 2 in the field
    FieldType num_0; // number 0 in the field
    int fieldByteSize;

    vector<byte> h;//a string accumulated that should be hashed in the comparing views function.



    vector<shared_ptr<ProtocolPartyData>> parties; // array of channels
    boost::asio::io_service io_service;
    ArithmeticCircuit circuit;
    vector<FieldType> gateShareArr; // my share of the gate (for all gates)
    vector<FieldType> alpha;
    vector<FieldType> alpha_for_mults;
    int index_for_mults = 0;
    vector<FieldType> random_for_inputs;

    vector<FieldType> r_for_verify_key1;
    vector<FieldType> r_for_comparing_views;
    vector<FieldType> random_element_for_verify;



    vector<FieldType> x_triple;
    vector<FieldType> y_triple;
    vector<FieldType> z_triple;
    vector <FieldType> beta_triple;

    int shareIndex; // number of shares for inputs
    int mult_count = 0;

    vector<int> myInputs;
    string s;

    shared_ptr<CommParty> leftChannel; // the channel with party i minus 1
    shared_ptr<CommParty> rightChannel; // the channel with party i plus 1

    TemplateField<FieldType>* field;

public:
    ProtocolParty(int argc, char* argv[]);

    /**
     * This method runs the protocol:
     * 2. Generate Randomness
     * 3. Input Preparation
     * 4. Computation Phase
     * 5. Verification Phase
     * 6. Output Phase
     */
    void run() override;

    virtual bool hasOffline() {
        return true;
    }


    virtual bool hasOnline() override {
        return true;
    }

    /**
     * This method runs the protocol:
     * Generate Randomness
     *
     */
    virtual void runOffline() override;

    /**
     * This method runs the protocol:
     * Input Preparation
     * Computation Phase
     * Verification Phase
     * Output Phase
     */
    virtual void runOnline() override;

    void roundFunctionSync(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs, int round);
    void exchangeData(vector<vector<byte>> &sendBufs,vector<vector<byte>> &recBufs, int first, int last);
    void sendNext(vector<byte> &sendBufs, vector<byte> &recBufs);

    /**
     * This method reads text file and inits a vector of Inputs according to the file.
     */
    void readMyInputs();

    /**
     * In case the user use the protocol with the offline/online mode, he should also update the iteration number.
     */
    void setIteration(int iteration) {
        this->iteration = iteration;
    }

    /**
     * We describe the protocol initialization.
     * In particular, some global variables are declared and initialized.
     */
    void initializationPhase();

    /**
     * this method prepares all the randomness which required in the protocol.
     */
    void generateRandomness31Bits();

    /**
     * this method prepares all the randomness which required in the protocol.
     */
    void generateRandomness61Bits();

    /**
     * this method prepares the c in each triple by multiplying a and b
     */
    void honestMultiplication(FieldType *a, FieldType *b, vector<FieldType> &cToFill, int numOfMult);

    /**
     * this method prepare the gate share array according the inputs and the shares.
     */
    void inputPreparation();

    /**
     * Walk through the circuit and evaluate the gates. Always take as many gates at once as possible,
     * i.e., all gates whose inputs are ready.
     * We first process all random gates, then alternately process addition and multiplication gates.
     */
    void computationPhase();

    /**
     * Process all additions which are ready.
     * Return number of processed gates.
     */
    int processNotMult();

    /**
     * Process all multiplications which are ready.
     * Return number of processed gates.
     */
    // int processMultiplications(vector<FieldType>& x, vector<FieldType>& y, vector<FieldType>& results, vector<FieldType>& alpha_for_mult, int begin_alpha,
    //                           int sizeOfSendBufsElements, int begin, int end);

    int processMultiplications();

    /**
     * this method check if the triples are currect.
     */
    void verification(int numOfMult);

    /**
     * this method open the shares and reconstructs the secrets.
     */
    void openShare(int numOfRandomShares, vector<FieldType> &shares, vector<FieldType> &secrets);

    /**
     * this method generate common key between the parties.
     */
    vector<byte> generateCommonKey(vector<FieldType>& random_elements);

    /**
     * this method sends the results from verification phase.
     */
    void comparingViews();

    /**
     * Walk through the circuit and reconstruct output gates.
     */
    void outputPhase();

    ~ProtocolParty();
};

template <class FieldType>
ProtocolParty<FieldType>::ProtocolParty(int argc, char* argv[]) : Protocol("ReplicatedSecretSharing3PartiesArithmetic", argc, argv) {
    CmdParser parser;
    this->times = stoi(parser.getValueByKey(arguments, "internalIterationsNumber"));
    vector<string> subTaskNames{"Offline", "preparationPhase", "Online", "inputPhase", "ComputePhase", "VerificationPhase", "outputPhase"};
    timer = new Measurement(*this, subTaskNames);

    string fieldType = parser.getValueByKey(arguments, "fieldType");
    if(fieldType.compare("ZpMersenne") == 0)
    {
        field = new TemplateField<FieldType>(2147483647);
    } else if(fieldType.compare("ZpMersenne61") == 0)
    {
        field = new TemplateField<FieldType>(0);
    }
    this->inputsFile = parser.getValueByKey(arguments, "inputFile");
    this->outputFile = parser.getValueByKey(arguments, "outputFile");

    m_partyId = stoi(parser.getValueByKey(arguments, "partyID"));
    s = to_string(m_partyId);
    circuit.readCircuit(parser.getValueByKey(arguments, "circuitFile").c_str());
    circuit.reArrangeCircuit();
    M = circuit.getNrOfGates();
    numOfInputGates = circuit.getNrOfInputGates();
    numOfOutputGates = circuit.getNrOfOutputGates();
    myInputs.resize(numOfInputGates);
    shareIndex = numOfInputGates;
    string partiesFilePath = parser.getValueByKey(arguments, "partiesFile");
    parties = MPCCommunication::setCommunication(io_service, m_partyId, N, partiesFilePath);

    int R = 0, L = 1; // TO DO: communication

    if(m_partyId == 1)
    {
        R = 1;
        L = 0;
    }

    leftChannel = parties[L]->getChannel();
    rightChannel = parties[R]->getChannel();

    string tmp = "init times";

    byte tmpBytes[20];
    for (int i=0; i<parties.size(); i++){
        if (parties[i]->getID() < m_partyId){
            parties[i]->getChannel()->write(tmp);
            parties[i]->getChannel()->read(tmpBytes, tmp.size());
        } else {
            parties[i]->getChannel()->read(tmpBytes, tmp.size());
            parties[i]->getChannel()->write(tmp);
        }
    }

    readMyInputs();
    initializationPhase();
}

template <class FieldType>
void ProtocolParty<FieldType>::readMyInputs()
{
    ifstream myfile;
    int input;
    int i =0;
    myfile.open(inputsFile);
    do {
        myfile >> input;
        myInputs[i] = input;
        i++;
    } while(!(myfile.eof()));
    myfile.close();

}

/**
     * Executes the protocol.
     */
template <class FieldType>
void ProtocolParty<FieldType>::run() {
    for (int i=0; i<times; i++){
        iteration = i;
        timer->startSubTask("Offline", iteration);
        runOffline();
        timer->endSubTask("Offline", iteration);
        timer->startSubTask("Online", iteration);
        runOnline();
        timer->endSubTask("Online", iteration);
    }
}

template <class FieldType>
void ProtocolParty<FieldType>::runOffline() {

    index_for_mults = 0;
    auto t1 = high_resolution_clock::now();
    timer->startSubTask("preparationPhase", iteration);

    if(fieldByteSize > 4) {
        generateRandomness61Bits();
    } else {
        generateRandomness31Bits();
    }

    timer->endSubTask("preparationPhase", iteration);
}

template <class FieldType>
void ProtocolParty<FieldType>::runOnline() {
    timer->startSubTask("inputPhase", iteration);
    inputPreparation();
    comparingViews();
    timer->endSubTask("inputPhase", iteration);

    timer->startSubTask("ComputePhase", iteration);
    computationPhase();
    timer->endSubTask("ComputePhase", iteration);

    timer->startSubTask("VerificationPhase", iteration);
    //calc the number of times we need to run the verification
    int iterations =   (5 + field->getElementSizeInBytes() - 1) / field->getElementSizeInBytes();
    if(circuit.getNrOfMultiplicationGates() > 0) {
        for(int i=0; i<iterations; i++) {
            verification(circuit.getNrOfMultiplicationGates() + circuit.getNrOfInputGates());
        }
    }
    timer->endSubTask("VerificationPhase", iteration);

    timer->startSubTask("outputPhase", iteration);
    outputPhase();
    timer->endSubTask("outputPhase", iteration);

}

template <class FieldType>
void ProtocolParty<FieldType>::initializationPhase()
{
    /**
     * The indexes in gate share arr:
     * for gates[k] :
     *      gateShareArr[2*k] = s
     *      gateShareArr[2*k+1] = t
     */
    gateShareArr.resize((M - circuit.getNrOfOutputGates())*4); // my share of the gate (for all gates)

    inv_3 = (field->GetElement(1))/(field->GetElement(3)); // calculate the inverse of 3 in the field

    fieldByteSize = field->getElementSizeInBytes();

    LEFT = (m_partyId - 1) % 3; // the number of party i Minus 1
    RIGHT = (m_partyId +1) % 3; // the number of party i Plus 1



    if(m_partyId == 0) {
        LEFT = 2;
        RIGHT = 1;

    }

    num_2 = (field->GetElement(2));

    num_0 = (field->GetElement(0));

    alpha.resize(2*circuit.getNrOfMultiplicationGates());
    random_for_inputs.resize(2 * circuit.getNrOfInputGates());

    if(fieldByteSize <= 4) {
        beta_triple.resize(2 * (circuit.getNrOfMultiplicationGates() + circuit.getNrOfInputGates()));
    } else {
        beta_triple.resize(circuit.getNrOfMultiplicationGates() + circuit.getNrOfInputGates());
    }

    r_for_verify_key1.resize(2* (16/field->getElementSizeInBytes() + 1)); // 16 bytes of aes key
    r_for_comparing_views.resize(2* (16/field->getElementSizeInBytes() + 1)); // 16 bytes of aes key
    random_element_for_verify.resize(2);//1 pair for the element

    x_triple.resize(2* circuit.getNrOfMultiplicationGates());
    y_triple.resize(2* circuit.getNrOfMultiplicationGates());
    z_triple.resize(2* circuit.getNrOfMultiplicationGates());

}

template <class FieldType>
void ProtocolParty<FieldType>::generateRandomness31Bits() {

    vector<byte> sendBufsBytes_1;
    int numOfMult = circuit.getNrOfMultiplicationGates();
    int rounds =  (16/field->getElementSizeInBytes() + 1);
    SecretKey key;
    vector<byte> key_i_1(16);

    int numOfTriples = numOfMult;

    if(fieldByteSize <= 4) {
        numOfTriples = numOfMult*2;
    }

    alpha_for_mults.resize(2*(numOfInputGates+1));

    PrgFromOpenSSLAES prg((numOfMult * 9 * fieldByteSize + (numOfInputGates+3*rounds + 3 )* fieldByteSize )/16 + 1);
    PrgFromOpenSSLAES prg2((numOfMult * 9 * fieldByteSize + (numOfInputGates+3*rounds + 3 )* fieldByteSize)/16 + 1);

    key = prg.generateKey(128);

    sendBufsBytes_1 = key.getEncoded();

    sendNext(sendBufsBytes_1,  key_i_1);

    SecretKey sk(key_i_1, "aes");

    prg2.setKey(sk); // ki-1

    prg.setKey(key); // ki

    for(int j=0; j<2*numOfMult;j++) {
        alpha[j] = field->GetElement(prg2.getRandom32()) - field->GetElement(prg.getRandom32());
    }

    for(int j=0; j<2*(numOfInputGates+1);j++) {
        alpha_for_mults[j] = field->GetElement(prg2.getRandom32()) - field->GetElement(prg.getRandom32());
    }

    // Generating Random Sharing For Inputs

    int numOfInputs = circuit.getNrOfInputGates();
    FieldType ri, riMinus1;

    for(int j=0; j < numOfInputs;j++) {
        ri = field->GetElement(prg.getRandom32());
        riMinus1 = field->GetElement(prg2.getRandom32());
        random_for_inputs[2*j] = riMinus1 - ri;
        random_for_inputs[2*j + 1] = num_0 - (num_2 * riMinus1) - ri;
    }

    // Generating Random Sharing For The Verification Stage and comparing views


    for(int j=0; j < rounds;j++) {
        ri = field->GetElement(prg.getRandom32());
        riMinus1 = field->GetElement(prg2.getRandom32());
        r_for_verify_key1[2*j] = riMinus1 - ri;
        r_for_verify_key1[2*j + 1] = num_0 - (num_2 * riMinus1) - ri;
    }


    for(int j=0; j < rounds;j++) {
        ri = field->GetElement(prg.getRandom32());
        riMinus1 = field->GetElement(prg2.getRandom32());
        r_for_comparing_views[2*j] = riMinus1 - ri;
        r_for_comparing_views[2*j + 1] = num_0 - (num_2 * riMinus1) - ri;
    }


    ri = field->GetElement(prg.getRandom32());
    riMinus1 = field->GetElement(prg2.getRandom32());
    random_element_for_verify[0] = riMinus1 - ri;
    random_element_for_verify[1] = num_0 - (num_2 * riMinus1) - ri;


}



template <class FieldType>
void ProtocolParty<FieldType>::generateRandomness61Bits() {

    vector<byte> sendBufsBytes_1;
    int numOfMult = circuit.getNrOfMultiplicationGates();
    SecretKey key;
    vector<byte> key_i_1(16);

    int numOfTriples = numOfMult;
    int rounds =  (16/field->getElementSizeInBytes() + 1);

    if(fieldByteSize <= 4) {
        numOfTriples = numOfMult*2;
    }

    alpha_for_mults.resize(numOfInputGates+1);

    PrgFromOpenSSLAES prg((numOfMult * 5 * fieldByteSize + (numOfInputGates+3*rounds + 2 )* fieldByteSize)/16 + 1);
    PrgFromOpenSSLAES prg2((numOfMult * 5 * fieldByteSize + (numOfInputGates+3*rounds + 2 ) * fieldByteSize)/16 + 1);

    key = prg.generateKey(128);

    sendBufsBytes_1 = key.getEncoded();

    sendNext(sendBufsBytes_1,  key_i_1);

    SecretKey sk(key_i_1, "aes");

    prg2.setKey(sk); // ki-1

    prg.setKey(key); // ki

    for(int j=0; j<2*numOfMult;j++) {
        alpha[j] = field->GetElement(prg2.getRandom64()) - field->GetElement(prg.getRandom64());
    }

    for(int j=0; j<numOfInputGates+1;j++) {
        alpha_for_mults[j] = field->GetElement(prg2.getRandom64()) - field->GetElement(prg.getRandom64());
    }

    // Generating Random Sharing For Inputs

    int numOfInputs = circuit.getNrOfInputGates();
    FieldType ri, riMinus1;

    for(int j=0; j < numOfInputs;j++) {
        ri = field->GetElement(prg.getRandom64());
        riMinus1 = field->GetElement(prg2.getRandom64());
        random_for_inputs[2*j] = riMinus1 - ri;
        random_for_inputs[2*j + 1] = num_0 - (num_2 * riMinus1) - ri;
    }

    // Generating Random Sharing For The Verification Stage


    for(int j=0; j < rounds;j++) {
        ri = field->GetElement(prg.getRandom64());
        riMinus1 = field->GetElement(prg2.getRandom64());
        r_for_verify_key1[2*j] = riMinus1 - ri;
        r_for_verify_key1[2*j + 1] = num_0 - (num_2 * riMinus1) - ri;
    }


    for(int j=0; j < rounds;j++) {
        ri = field->GetElement(prg.getRandom64());
        riMinus1 = field->GetElement(prg2.getRandom64());
        r_for_comparing_views[2*j] = riMinus1 - ri;
        r_for_comparing_views[2*j + 1] = num_0 - (num_2 * riMinus1) - ri;
    }

    ri = field->GetElement(prg.getRandom64());
    riMinus1 = field->GetElement(prg2.getRandom64());
    random_element_for_verify[0] = riMinus1 - ri;
    random_element_for_verify[1] = num_0 - (num_2 * riMinus1) - ri;


}


template <class FieldType>
void ProtocolParty<FieldType>::honestMultiplication(FieldType *a, FieldType *b, vector<FieldType> &cToFill, int numOfMult) {

    FieldType p2, d2;
    FieldType ri, riMinus1;
    vector<FieldType> sendBufsElements(numOfMult);
    vector<byte> sendBufsBytes(numOfMult*field->getElementSizeInBytes());
    vector<byte> recBufsBytes(numOfMult*field->getElementSizeInBytes());

    for(int k = 0; k < numOfMult ; k++)
    {

        ri = (a[(k * 2) + 1] * b[(k * 2) + 1] -
              (a[(k * 2)] * b[(k * 2)]) + alpha_for_mults[index_for_mults + k]) * inv_3;

        //send ri to pi+1 = RIGHT
        sendBufsElements[k] = ri;


    }

    //convert to bytes

//    for(int j=0; j < numOfMult;j++) {
//        field->elementToBytes(sendBufsBytes.data() + (j * fieldByteSize), sendBufsElements[j]);
//    }

    field->elementVectorToByteVector(sendBufsElements, sendBufsBytes);

    sendNext(sendBufsBytes, recBufsBytes);

    int fieldBytesSize = field->getElementSizeInBytes();

    for(int k = 0; k < numOfMult; k++)
    {
        riMinus1 = field->bytesToElement(recBufsBytes.data() + (k * fieldBytesSize));

        ri = sendBufsElements[k];

        cToFill[k * 2] = riMinus1 - ri; // ei

        cToFill[k * 2 + 1] = num_0 - (num_2 * riMinus1) - ri; // fi
    }

    index_for_mults+=numOfMult;
}




template <class FieldType>
void ProtocolParty<FieldType>::verification(int numOfMult)
{
    PrgFromOpenSSLAES prg(numOfMult/4+1);
   vector<byte> keyVector(16);

    keyVector = generateCommonKey(r_for_verify_key1);
    // generating 128 bit AES key
    SecretKey sk(keyVector, "aes");
    prg.setKey(sk);


    if(fieldByteSize > 4) {
        for(int j=0; j<numOfMult;j++) {

            beta_triple[j] = field->GetElement((prg.getRandom64() >> 3));
        }
    } else {
        for(int j=0; j<numOfMult;j++) {

            beta_triple[j] = field->GetElement(prg.getRandom32());
        }
    }


    //preapre x,y,z for the verification sub protocol
    vector<FieldType> neededShares(numOfMult*4);


    int index = 0;
    for (int k = 0; k < numOfInputGates; k++) {

        auto gate = circuit.getGates()[k];

        if (gate.gateType == INPUT) {
            neededShares[4*index] = gateShareArr[gate.output*4];
            neededShares[4*index+1] = gateShareArr[gate.output*4+1];
            neededShares[4*index+2] = gateShareArr[gate.output*4+2];
            neededShares[4*index+3] = gateShareArr[gate.output*4+3];

            index++;
        }
    }

    for (int k = numOfInputGates - 1; k < M - numOfOutputGates + 1; k++) {

        auto gate = circuit.getGates()[k];
        if (gate.gateType == MULT) {
            neededShares[4*index] = gateShareArr[gate.output*4];
            neededShares[4*index+1] = gateShareArr[gate.output*4+1];
            neededShares[4*index+2] = gateShareArr[gate.output*4+2];
            neededShares[4*index+3] = gateShareArr[gate.output*4+3];
            index++;
        }


    }
    vector<FieldType> u(2);
    vector<FieldType> w(2);
    vector<FieldType> ru(2);
    vector<FieldType> T(2);


    for(int i=0;i<numOfMult;i++){

        u[0] += beta_triple[i]*neededShares[i*4];
        u[1] += beta_triple[i]*neededShares[i*4 + 1];
        w[0] += beta_triple[i]*neededShares[i*4+2];
        w[1] += beta_triple[i]*neededShares[i*4+3];
    }


    //run the semi honest multiplication on u and r to get ru
    honestMultiplication(u.data(), random_element_for_verify.data(),ru, 1);

    T[0] = w[0] - ru[0];
    T[1] = w[1] - ru[1];


    //open [T]

    vector<FieldType> secretArr(1);

    openShare(1,T,secretArr);

    //check that T=0
    if(secretArr[0] != *field->GetZero()) {
        if(flag_print)
            cout<<"bassssssaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"<<endl;
        return;
    }
    else {
        if(flag_print)
            cout<<"yessssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss"<<endl;
        return;
    }

}

template <class FieldType>
void ProtocolParty<FieldType>::comparingViews()
{

    vector<byte> keyVector;
    unsigned int hashSize = 16;
    unsigned char digest[hashSize];
    vector<byte> recBufBytes(16);
    vector<byte> sendBuf;

    unsigned char hashedViewsForCompare[hashSize];

    //not sure this is even needed, why should the aes key be generated now and not before the entire computation
    keyVector = generateCommonKey(r_for_comparing_views);

    // gcm initialization vector
    unsigned char iv1[] = {0xe0, 0xe0, 0x0f, 0x19, 0xfe, 0xd7, 0xba, 0x01,
                           0x36, 0xa7, 0x97, 0xf3};

    unsigned char *key = reinterpret_cast<unsigned char*>(keyVector.data());


    HashEncrypt hashObj = HashEncrypt(key, iv1, 12);

    hashObj.getHashedDataOnce(reinterpret_cast<unsigned char*> (h.data()), h.size(), digest, &hashSize);



    //copy the digest to each send buf
    copy_byte_array_to_byte_vector(digest,hashSize,sendBuf,0);
    sendNext(sendBuf, recBufBytes);

    // check the result

    bool flag = true;

    for(int j = 0; j < 16;j++) {
        if (recBufBytes[j] !=  digest[j]) {
            cout << "cheating hi" << endl;
            flag = false;
        }
    }


    if(flag) {
        cout << "no cheating" << endl;
    }

}

template <class FieldType>
void ProtocolParty<FieldType>::inputPreparation()
{
    int robin = 0;
    // the number of random double sharings we need altogether
    vector<FieldType> x1(N),y1(N);
    vector<vector<FieldType>> sendBufsElements(N);
    vector<vector<byte>> sendBufsBytes(N);
    vector<vector<byte>> recBufBytes(N);
    vector<vector<FieldType>> recBufElements(N);
    vector<FieldType> inputShares(circuit.getNrOfInputGates()*2);

    int input;
    int index = 0;
    vector<int> sizes(N);
    FieldType s1,s2,s3,t1,t2,t3, r1,r2,r3;

    // prepare the shares for the inputs
    for (int k = 0; k < numOfInputGates; k++)
    {
        if(circuit.getGates()[k].gateType == INPUT) {
            //get the expected sized from the other parties
            sizes[(circuit.getGates()[k].party)]++;  // MEITAL

            // send to party (which need this gate) your share for this gate
            sendBufsElements[circuit.getGates()[k].party].push_back(random_for_inputs[2*k]);// send t
        }
    }

    for(int i=0; i < N; i++)
    {

        sendBufsBytes[i].resize(sendBufsElements[i].size()*fieldByteSize);
        recBufBytes[i].resize(sizes[m_partyId]*fieldByteSize);
//        for(int j=0; j<sendBufsElements[i].size();j++) {
//            field->elementToBytes(sendBufsBytes[i].data() + (j * fieldByteSize), sendBufsElements[i][j]);
//        }

        field->elementVectorToByteVector(sendBufsElements[i], sendBufsBytes[i]);
    }

    roundFunctionSync(sendBufsBytes, recBufBytes,10);


    //turn the bytes to elements
    for(int i=0; i < N; i++)
    {
        recBufElements[i].resize(((recBufBytes[i].size()) / fieldByteSize));
        for(int j=0; j<recBufElements[i].size();j++) {
            recBufElements[i][j] = field->bytesToElement(recBufBytes[i].data() + ( j * fieldByteSize));
        }
    }

    vector<int> counters(N);

    for(int i=0; i<N; i++){
        counters[i] =0;
    }

    int counter = 0;
    FieldType t_i_1, ti, m_t, m_s;

    vector<vector<FieldType>> sendBufsElements_1(N);
    vector<vector<byte>> sendBufsBytes_1(N);
    vector<vector<byte>> recBufBytes_1(N);
    vector<vector<FieldType>> recBufElements_1(N);

    index = 0;
    vector<int> sizes_1(N);
    for(int i=0; i<N; i++){
        sizes_1[i] = 0;
    }

     //  The parties run reconstruct
    for (int k = 0; k < numOfInputGates; k++)
    {
        if(circuit.getGates()[k].gateType == INPUT && circuit.getGates()[k].party == m_partyId)
        {
            t_i_1 = field->bytesToElement(recBufBytes[LEFT].data() + (counter*fieldByteSize));
            ti = field->bytesToElement(recBufBytes[RIGHT].data() + ((counter)*fieldByteSize));

            m_t = random_for_inputs[k*2];
            m_s = random_for_inputs[k*2 + 1];
            FieldType value = t_i_1 + ti + m_t;
            // my output: reconstruct received shares
            if (value != *field->GetZero())
            {
                // someone cheated!
                cout << "cheating!!!" << '\n';

                return;
            }

            FieldType v = t_i_1 - m_s; // holding r

            auto input = myInputs[index];
            index++;

            sendBufsElements_1[LEFT].push_back(field->GetElement(input) - v);
            sendBufsElements_1[RIGHT].push_back(field->GetElement(input) - v);
            sendBufsElements_1[m_partyId].push_back(field->GetElement(input) - v);

            counter++;

            gateShareArr[4*(circuit.getGates()[k].output)] = m_t;// set the share sent from the party owning the input
            inputShares[2*k] = m_t;
            gateShareArr[4*(circuit.getGates()[k].output)+1] = m_s - (field->GetElement(input) - v);
            inputShares[2*k+1] = gateShareArr[4*(circuit.getGates()[k].output)+1];

        }
        sizes_1[circuit.getGates()[k].party]++;
    }

    for(int i=0; i < N; i++)
    {
        sendBufsBytes_1[i].resize(sendBufsElements_1[i].size()*fieldByteSize);
        recBufBytes_1[i].resize(sizes_1[i]*fieldByteSize); // MEITAL!!
//        for(int j=0; j<sendBufsElements_1[i].size();j++) {
//            field->elementToBytes(sendBufsBytes_1[i].data() + (j * fieldByteSize), sendBufsElements_1[i][j]);
//        }

        field->elementVectorToByteVector(sendBufsElements_1[i], sendBufsBytes_1[i]);
    }

    roundFunctionSync(sendBufsBytes_1, recBufBytes_1,11);

    //turn the bytes to elements
    for(int i=0; i < N; i++)
    {
        recBufElements_1[i].resize(((recBufBytes_1[i].size()) / fieldByteSize));
        for(int j=0; j<recBufElements_1[i].size();j++) {
            recBufElements_1[i][j] = field->bytesToElement(recBufBytes_1[i].data() + ( j * fieldByteSize));
        }
    }

    //update the empty string to include the b's recieved from the other parties.

    //move the first set of b's of party 0, the other 2 vectors will need to be copied
    h = move(recBufBytes_1[0]);

    for(int i=1; i<N; i++)
    {
        h.insert(h.end(), recBufBytes_1[i].begin(), recBufBytes_1[i].end());
    }

    vector<int> counters_1(N);

    for(int i=0; i<N; i++){
        counters_1[i] = 0;
    }

    //every party compute locally the shares of inputs

    for (int k = 0; k < numOfInputGates; k++) {
        if (circuit.getGates()[k].gateType == INPUT && circuit.getGates()[k].party != m_partyId) {
            gateShareArr[4*(circuit.getGates()[k].output)] = random_for_inputs[2*k];// set the share sent from the party owning the input
            inputShares[2*k] = gateShareArr[4*(circuit.getGates()[k].output)];
            gateShareArr[4*(circuit.getGates()[k].output)+1] = random_for_inputs[2*k + 1] -
                    recBufElements_1[circuit.getGates()[k].party][counters_1[circuit.getGates()[k].party]];
            inputShares[2*k+1] = gateShareArr[4*(circuit.getGates()[k].output)+1];
            counters_1[circuit.getGates()[k].party] += 1;


        }
    }
    //set this random share to an entire array so we can use the semi honest multiplication
    vector<FieldType> resultMult(numOfInputGates*2);

    for(int i=0; i<numOfInputGates; i++){
        resultMult[2*i] = random_element_for_verify[0];
        resultMult[2*i+1] = random_element_for_verify[1];

    }

    //run the semi honest multiplication to get the second part of each share
    honestMultiplication(inputShares.data(), resultMult.data(),resultMult, numOfInputGates);

    //set the resulted multiplication to the array of shares

    for (int k = 0; k < numOfInputGates; k++)
    {
        if(circuit.getGates()[k].gateType == INPUT)
        {
            //set the second part of the share
            gateShareArr[circuit.getGates()[k].output*4+2] = resultMult[2*k];
            gateShareArr[circuit.getGates()[k].output*4+3] = resultMult[2*k+1];

        }
    }

}

template <class FieldType>
void ProtocolParty<FieldType>::openShare(int numOfRandomShares, vector<FieldType> &shares, vector<FieldType> &secrets){

    int robin = 0;

    vector<vector<FieldType>> sendBufsElements(N);
    vector<vector<byte>> sendBufsBytes(N);
    vector<vector<byte>> recBufBytes(N);
    vector<vector<FieldType>> recBufElements(N);

    int input;
    int index = 0;
    FieldType s1,s2,s3,t1,t2,t3, r1,r2,r3;

    for (int k = 0; k < numOfRandomShares; k++)
    {
        sendBufsElements[RIGHT].push_back(shares[2*k]);// t
    }

    int fieldByteSize = field->getElementSizeInBytes();

    for(int i=0; i < N; i++)
    {
        sendBufsBytes[i].resize(sendBufsElements[i].size()*fieldByteSize);

//        for(int j=0; j<sendBufsElements[i].size();j++) {
//            field->elementToBytes(sendBufsBytes[i].data() + (j * fieldByteSize), sendBufsElements[i][j]);
//        }

        field->elementVectorToByteVector(sendBufsElements[i], sendBufsBytes[i]);
    }

    recBufBytes[LEFT].resize(numOfRandomShares*fieldByteSize);
    recBufBytes[RIGHT].resize(0);
    recBufBytes[m_partyId].resize(numOfRandomShares*fieldByteSize);

    roundFunctionSync(sendBufsBytes, recBufBytes,10);


    //turn the bytes to elements
    for(int i=0; i < N; i++)
    {
        recBufElements[i].resize(((recBufBytes[i].size()) / fieldByteSize));
        for(int j=0; j<recBufElements[i].size();j++) {
            recBufElements[i][j] = field->bytesToElement(recBufBytes[i].data() + ( j * fieldByteSize));
        }
    }

    int counter = 0;
    FieldType t_i_1, ti, m_t, m_s;

    vector<vector<FieldType>> sendBufsElements_1(N);
    vector<vector<byte>> sendBufsBytes_1(N);
    vector<vector<byte>> recBufBytes_1(N);
    vector<vector<FieldType>> recBufElements_1(N);

    index = 0;
    vector<int> sizes_1(N);
    for(int i=0; i<N; i++){
        sizes_1[i] = 0;
    }

    //  The parties run reconstruct

    for (int k = 0; k < numOfRandomShares; k++) {
        t_i_1 = field->bytesToElement(recBufBytes[LEFT].data() + (counter * fieldByteSize));

        m_s = shares[k * 2 + 1];

        FieldType v = t_i_1 - m_s; // holding r

        secrets[k] = v;

        counter++;
    }
}

template <class FieldType>
vector<byte> ProtocolParty<FieldType>::generateCommonKey(vector<FieldType>& random_elements){

    //calc the number of elements needed for 128 bit AES key
   int numOfRandomShares = 16/field->getElementSizeInBytes() + 1;
    vector<FieldType> aesArray(numOfRandomShares);

   // vector<FieldType> aesArray(numOfRandomShares);
    vector<byte> aesKey(numOfRandomShares*fieldByteSize);

    //generate enough random shares for the AES key
    openShare(numOfRandomShares, random_elements, aesArray);

    //turn the aes array into bytes to get the common aes key.
    for(int i=0; i<numOfRandomShares;i++){

        for(int j=0; j<numOfRandomShares;j++) {
            field->elementToBytes(aesKey.data() + (j * fieldByteSize), aesArray[j]);
        }
    }

    //reduce the size of the key to 16 bytes
    aesKey.resize(16);

    return aesKey;

}


template <class FieldType>
void ProtocolParty<FieldType>::computationPhase() {

    currentCirciutLayer = 0;
    mult_count = 0;

    int count = 0, num;
    int c_currentCirciutLayer, c_nextCirciutLayer;

    int numOfLayers = circuit.getLayers().size();
    for(int i=0; i<numOfLayers-1;i++){

        currentCirciutLayer = i;
        count = processNotMult();

        c_currentCirciutLayer = circuit.getLayers()[currentCirciutLayer];
        c_nextCirciutLayer =  circuit.getLayers()[currentCirciutLayer+1];

        count += processMultiplications();
    }
}


template <class FieldType>
int ProtocolParty<FieldType>::processNotMult(){
    int count=0;
    for(int k=circuit.getLayers()[currentCirciutLayer]; k < circuit.getLayers()[currentCirciutLayer+1]; k++)
    {
        // add gate
        if(circuit.getGates()[k].gateType == ADD)
        {
            gateShareArr[circuit.getGates()[k].output * 4] = gateShareArr[circuit.getGates()[k].input1 * 4] + gateShareArr[circuit.getGates()[k].input2 * 4]; // t+t
            gateShareArr[(circuit.getGates()[k].output * 4) + 1] = gateShareArr[(circuit.getGates()[k].input1 * 4) + 1] + gateShareArr[(circuit.getGates()[k].input2 * 4) + 1]; // s+s

            gateShareArr[(circuit.getGates()[k].output * 4) + 2] = gateShareArr[(circuit.getGates()[k].input1 * 4) + 2] + gateShareArr[(circuit.getGates()[k].input2 * 4) + 2]; // t+t
            gateShareArr[(circuit.getGates()[k].output * 4) + 3] = gateShareArr[(circuit.getGates()[k].input1 * 4) + 3] + gateShareArr[(circuit.getGates()[k].input2 * 4) + 3]; // s+s

            count++;

        }
        else if(circuit.getGates()[k].gateType == SCALAR)
        {
            long scalar(circuit.getGates()[k].input2);
            FieldType e = field->GetElement(scalar);

            gateShareArr[circuit.getGates()[k].output * 4] = gateShareArr[circuit.getGates()[k].input1 * 4]; // t
            gateShareArr[(circuit.getGates()[k].output * 4) + 1] = gateShareArr[(circuit.getGates()[k].input1 * 4) + 1] * e; // s*e
            gateShareArr[(circuit.getGates()[k].output * 4) + 2] = gateShareArr[(circuit.getGates()[k].input1 * 4) + 2] * e; // s*e
            gateShareArr[(circuit.getGates()[k].output * 4) + 3] = gateShareArr[(circuit.getGates()[k].input1 * 4) + 3] * e; // s*e

            count++;
        }

        else if(circuit.getGates()[k].gateType == SCALAR_ADD)
        {
            long scalar(circuit.getGates()[k].input2);
            FieldType e = field->GetElement(scalar);

            gateShareArr[circuit.getGates()[k].output * 4] = gateShareArr[circuit.getGates()[k].input1 * 4]; // t
            gateShareArr[(circuit.getGates()[k].output * 4) + 1] = gateShareArr[(circuit.getGates()[k].input1 * 4) + 1] - e; // s-e
            gateShareArr[(circuit.getGates()[k].output * 4) + 2] = gateShareArr[(circuit.getGates()[k].input1 * 4) + 2] ; // t
            gateShareArr[(circuit.getGates()[k].output * 4) + 3] = gateShareArr[(circuit.getGates()[k].input1 * 4) + 3] - e; // s-e

            count++;
        }
    }

    return count;

}

/**
 * the Function process all multiplications which are ready.
 * @return the number of processed gates.
 */
template <class FieldType>
int ProtocolParty<FieldType>::processMultiplications()
{
    int last = circuit.getLayers()[currentCirciutLayer+1];
    int first = circuit.getLayers()[currentCirciutLayer];
    int size = (last - first);
    int index = 0;
    FieldType p2, d2;
    FieldType ri, r_i_1;
    vector<FieldType> sendBufsElements(2*size);
    vector<byte> sendBufsBytes(2*size*fieldByteSize);
    vector<byte> recBufsBytes(2*size*fieldByteSize);

    for(int k = circuit.getLayers()[currentCirciutLayer]; k < circuit.getLayers()[currentCirciutLayer+1] ; k++)//go over only the logit gates
    {
        // its a multiplication which not yet processed and ready
        if(circuit.getGates()[k].gateType == MULT)
        {

            ri = (gateShareArr[(circuit.getGates()[k].input1 * 4) + 1] * gateShareArr[(circuit.getGates()[k].input2 * 4) + 1] -
                    gateShareArr[(circuit.getGates()[k].input1 * 4)] * gateShareArr[(circuit.getGates()[k].input2 * 4)] + alpha[mult_count]) * inv_3;

            //send ri to pi+1 = RIGHT
            sendBufsElements[2*index] = ri;
            mult_count++;

            ri = (gateShareArr[(circuit.getGates()[k].input1 * 4) + 3] * gateShareArr[(circuit.getGates()[k].input2 * 4) + 1] -
                  gateShareArr[(circuit.getGates()[k].input1 * 4+2)] * gateShareArr[(circuit.getGates()[k].input2 * 4)] + alpha[mult_count]) * inv_3;

            //send ri to pi+1 = RIGHT
            sendBufsElements[2*index + 1] = ri;

            index++;
            mult_count++;
        }
    }

    //convert to bytes
//    for(int j=0; j < size; j++) {
//        field->elementToBytes(sendBufsBytes.data() + (j * fieldByteSize), sendBufsElements[j]);
//    }

    field->elementVectorToByteVector(sendBufsElements, sendBufsBytes);

    sendNext(sendBufsBytes, recBufsBytes);

    index = 0;

    for(int k = first; k < last; k++) {

        if(circuit.getGates()[k].gateType == MULT) {

                r_i_1 = field->bytesToElement(recBufsBytes.data() + (2*index * fieldByteSize));

                ri = sendBufsElements[2*index];

                gateShareArr[(circuit.getGates()[k].output) * 4] = r_i_1 - ri; // ei

                gateShareArr[(circuit.getGates()[k].output) * 4 + 1] = num_0 - (num_2 * r_i_1) - ri; //  fi

                r_i_1 = field->bytesToElement(recBufsBytes.data() + ((2*index +1) * fieldByteSize));

                ri = sendBufsElements[2*index+1];

                gateShareArr[(circuit.getGates()[k].output) * 4+2] = r_i_1 - ri; // ei

                gateShareArr[(circuit.getGates()[k].output) * 4 + 3] = num_0 - (num_2 * r_i_1) - ri; //  fi


            index++;
        }

    }


    return index;
}

/**
 * the function Walk through the circuit and reconstruct output gates.
 * @param circuit
 * @param gateShareArr
 * @param alpha
 */
template <class FieldType>
void ProtocolParty<FieldType>::outputPhase()
{
    int count=0;
    vector<FieldType> x1(N); // vector for the shares of my outputs
    vector<vector<FieldType>> sendBufsElements(N);
    vector<vector<byte>> sendBufsBytes(N);
    vector<vector<byte>> recBufBytes(N);

    FieldType num;
    ofstream myfile;
    myfile.open(outputFile);

    for(int k=M-numOfOutputGates; k < M; k++)
    {
        if(circuit.getGates()[k].gateType == OUTPUT)
        {
            // send to party (which need this gate) your share for this gate
            sendBufsElements[circuit.getGates()[k].party].push_back(gateShareArr[4*circuit.getGates()[k].input1]);// send t
        }
    }

    for(int i=0; i < N; i++)
    {
        sendBufsBytes[i].resize(sendBufsElements[i].size()*fieldByteSize);
        recBufBytes[i].resize(sendBufsElements[m_partyId].size()*fieldByteSize);
//        for(int j=0; j<sendBufsElements[i].size();j++) {
//            field->elementToBytes(sendBufsBytes[i].data() + (j * fieldByteSize), sendBufsElements[i][j]);
//        }

        field->elementVectorToByteVector(sendBufsElements[i], sendBufsBytes[i]);
    }

    roundFunctionSync(sendBufsBytes, recBufBytes,7);

    FieldType t_i_1, ti, m_t, m_s;
    int counter = 0;

    for(int k=M-numOfOutputGates ; k < M; k++) {
        if(circuit.getGates()[k].gateType == OUTPUT && circuit.getGates()[k].party == m_partyId)
        {

            t_i_1 = field->bytesToElement(recBufBytes[LEFT].data() + (counter*fieldByteSize));
            ti = field->bytesToElement(recBufBytes[RIGHT].data() + ((counter)*fieldByteSize));

            m_t = gateShareArr[(circuit.getGates()[k].input1)*4];
            m_s = gateShareArr[(circuit.getGates()[k].input1)*4 + 1];
            FieldType value = t_i_1 + ti + m_t;

            // my output: reconstruct received shares
            if (value != *field->GetZero())
            {
                // someone cheated!
                cout << "cheating!!!" << '\n';

                return;
            }

            FieldType v = t_i_1 - m_s;
            cout << "the result for "<< circuit.getGates()[k].input1 << " is : " << field->elementToString(t_i_1 - m_s) << '\n';

            counter++;
        }
    }

    // close output file
    myfile.close();
}

/**
 * communication
 * @tparam FieldType
 * @param sendBufs
 * @param recBufs
 * @param round
 */
template <class FieldType>
void ProtocolParty<FieldType>::roundFunctionSync(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs, int round) {

    int numThreads = 2;
    int numPartiesForEachThread = 1;

    recBufs[m_partyId] = sendBufs[m_partyId];
    //recieve the data using threads
    vector<thread> threads(numThreads);
    for (int t=0; t<numThreads; t++) {
        if ((t + 1) * numPartiesForEachThread <= parties.size()) {
            threads[t] = thread(&ProtocolParty::exchangeData, this, ref(sendBufs), ref(recBufs),
                                t * numPartiesForEachThread, (t + 1) * numPartiesForEachThread);
        } else {
            threads[t] = thread(&ProtocolParty::exchangeData, this, ref(sendBufs), ref(recBufs), t * numPartiesForEachThread, parties.size());
        }
    }
    for (int t=0; t<numThreads; t++){
        threads[t].join();
    }

}

template <class FieldType>
void ProtocolParty<FieldType>::sendNext(vector<byte> &sendBufs, vector<byte> &recBufs) {

    if (m_partyId == 0)
    {
        rightChannel->write(sendBufs.data(), sendBufs.size()); // write to party 2

        leftChannel->read(recBufs.data(), recBufs.size()); // read from party 3

    } else if (m_partyId == 1)
    {
        leftChannel->read(recBufs.data(), recBufs.size()); // read from party 1

        rightChannel->write(sendBufs.data(), sendBufs.size()); // write to party 3

    } else {

        rightChannel->write(sendBufs.data(), sendBufs.size()); // write to party 1

        leftChannel->read(recBufs.data(), recBufs.size()); // read from party 2
    }

}

template <class FieldType>
void ProtocolParty<FieldType>::exchangeData(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs, int first, int last){

    //cout<<"in exchangeData";
    for (int i=first; i < last; i++) {

        if ((m_partyId) < parties[i]->getID()) {


            //send shares to my input bits
            parties[i]->getChannel()->write(sendBufs[parties[i]->getID()].data(), sendBufs[parties[i]->getID()].size());
            //cout<<"write the data:: my Id = " << m_partyId - 1<< "other ID = "<< parties[i]->getID() <<endl;


            //receive shares from the other party and set them in the shares array
            parties[i]->getChannel()->read(recBufs[parties[i]->getID()].data(), recBufs[parties[i]->getID()].size());
            //cout<<"read the data:: my Id = " << m_partyId-1<< "other ID = "<< parties[i]->getID()<<endl;

        } else{


            //receive shares from the other party and set them in the shares array
            parties[i]->getChannel()->read(recBufs[parties[i]->getID()].data(), recBufs[parties[i]->getID()].size());
            //cout<<"read the data:: my Id = " << m_partyId-1<< "other ID = "<< parties[i]->getID()<<endl;



            //send shares to my input bits
            parties[i]->getChannel()->write(sendBufs[parties[i]->getID()].data(), sendBufs[parties[i]->getID()].size());
            //cout<<"write the data:: my Id = " << m_partyId-1<< "other ID = "<< parties[i]->getID() <<endl;


        }

    }

}


template <class FieldType>
ProtocolParty<FieldType>::~ProtocolParty()
{
    delete timer;
}

#endif /* PROTOCOL_H_ */
