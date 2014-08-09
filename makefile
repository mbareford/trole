SRCS = ${shell ls ./src/*.cpp}
OBJS = ${SRCS:./src/%.cpp=./obj/%.o}
LIBS = -lgsl -lgslcblas

CC = g++
#debug
CFLAGS = -O0 -g3 -Wall -fmessage-length=0
#release
#CFLAGS = -O1 -Wall -fmessage-length=0


./bin/trole: $(OBJS)
	$(CC) -o "./bin/trole" $(OBJS) $(LIBS)

./obj/%.o: ./src/%.cpp
	$(CC) $(CFLAGS) -c -o "$@" "$<"


clean:
	rm -f $(OBJS)
	rm -f ./bin/trole

info:
	@echo $(OBJS)
	@echo $(SRCS)
