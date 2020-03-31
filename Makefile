EXECUTABLE=kry.exe
SOURCE=kry.cpp 

compile:
	g++ -std=c++11  ${SOURCE} -o ${EXECUTABLE} -lgmp -lgmpxx

gen: compile
	./${EXECUTABLE} -g 64

enc: compile
	./${EXECUTABLE} -e 0xd 0x4d 0x2a

dec: compile
	./${EXECUTABLE} -d 0x25 0x4d 0xe

bre: compile
	./${EXECUTABLE} -b 0x3 0x2f72f3e766a6ee4497ec836160e05363be3c0eb5adbdddb47 0x850e13a7c8826f9a8682f611dffe8cadd75e57a351ab9864