cmake_minimum_required(VERSION 2.8)
project(ReplicatedSecretSharing3PartiesArithmeticNoTriples)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -maes -mpclmul -mbmi2")


INCLUDE_DIRECTORIES($ENV{HOME}/ntl-9.10.0/include $ENV{HOME} $ENV{HOME}/libscapi/install/include
      $ENV{HOME}/libscapi/lib/OTExtensionBristol )

set(SOURCE_FILES ProtocolParty.h main.cpp HashEncrypt.h HashEncrypt.cpp ProtocolTimer.h ProtocolTimer.cpp )


link_directories($ENV{HOME}/libscapi/install/lib)


add_executable(ReplicatedSecretSharing3PartiesArithmeticNoTriples ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(ReplicatedSecretSharing3PartiesArithmeticNoTriples $ENV{HOME}/ntl-9.10.0/src/ntl.a
		$ENV{HOME}/libscapi/libscapi.a gmp
		OTExtensionBristol libsimpleot.a libjsoncpp.a
        libboost_system.a libboost_thread.a libssl.a libcrypto.a pthread dl  z)
