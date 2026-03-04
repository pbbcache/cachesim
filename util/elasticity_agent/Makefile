CC = gcc
CFLAGS = -Wall -Wextra -g
LDFLAGS = 

SRCS = elasticity_agent.c
OBJS = $(SRCS:.c=.o)
DEPS = schedctl.h

TARGET = elasticity_agent

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS)
	cp $(TARGET) ../../simulator/test

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJS) $(TARGET)
