INCLDIRS = ./
CC = gcc
SRC += $(shell ls *.c)
OBJ += $(SRC:.c=.o)
CFLAGS += -O2 -Wall -g -c
LDFLAGS +=-lm
CPPFLAGS += -I$(INCLDIRS)
OUTPUT = $(shell pwd)/genetic-algorithms-run

$(OUTPUT): $(OBJ)
	$(CC) -o $@ $(OBJ) $(LDFLAGS)

-include $(SRC:%.c=%.P)

%.o : %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -MD -o $@ $<
	@cp $*.d $*.P
	@sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	     -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P
	@rm -f $*.d

.PHONY: clean
clean:
	rm -rf $(OUTPUT) $(OBJ) *.P
