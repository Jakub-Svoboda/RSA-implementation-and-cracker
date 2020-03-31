EXECUTABLE=kry.exe
SOURCE=kry.cpp 

compile:
	g++ -std=c++11  ${SOURCE} -o ${EXECUTABLE} -lgmp -lgmpxx

gen: compile
	./${EXECUTABLE} -g 8

enc: compile
	./${EXECUTABLE} -e 0x66fd30 0x66fd50 0x42

dec: compile
	./${EXECUTABLE} -d 0x66fd20 0x66fd50 0x66fd70

bre: compile
	./${EXECUTABLE} -b 0x3 0x2f72f3e766a6ee4497ec836160e05363be3c0eb5adbdddb47 0x850e13a7c8826f9a8682f611dffe8cadd75e57a351ab9864