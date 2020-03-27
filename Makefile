EXECUTABLE=kry.exe
SOURCE=kry.cpp 

compile:
	g++ ${SOURCE} -o ${EXECUTABLE}

gen: compile
	./${EXECUTABLE} -g B

enc: compile
	./${EXECUTABLE} -e E N M

dec: compile
	./${EXECUTABLE} -d D N C

bre: compile
	./${EXECUTABLE} -b E N C