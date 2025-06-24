##############################################################################
#   PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
#   Copyright (C) 2007-2016 University of Tennessee / UT-Battelle
#
#   Code Author:  Doug Hyatt
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
##############################################################################

SHELL   = /bin/sh
CC      = gcc

CFLAGS  += -pedantic -Wall -O3 -DSUPPORT_GZIP_COMPRESSED
LFLAGS = -lm -Wl,-O2 $(LDFLAGS) -lz

TARGET  = prodigal
TABLEUTIL = prodigal-table
ZTARGET  = zprodigal
SOURCES = $(shell echo *.c)
HEADERS = $(shell echo *.h) table_baked.h
OBJECTS = $(SOURCES:.c=.o)
ZOBJECTS = $(SOURCES:.c=.oz)

INSTALLDIR  = /usr/local/bin

all: $(TARGET) $(TABLEUTIL)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

$(TABLEUTIL): table.c table.h table_baked.o
	$(CC) $(CFLAGS) -DTABLE_UTIL -Wno-parentheses -o $@ $< table_baked.o $(LFLAGS)

table_baked.c table_baked.h: gc.prt maketable.sh
	./maketable.sh < gc.prt

# Too big to recompile on every header change
training_baked.o: training_baked.c training.h training_baked.h
	$(CC) $(CFLAGS) -c -o $@ $<

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o $@ $<

install: $(TARGET) $(TABLEUTIL)
	install -d -m 0755 $(INSTALLDIR)
	install -m 0755 $(TARGET) $(INSTALLDIR)
	install -m 0755 $(TABLEUTIL) $(INSTALLDIR)
 
uninstall:
	-rm $(INSTALLDIR)/$(TARGET)
	-rm $(INSTALLDIR)/$(TABLEUTIL)

clean:
	-rm -f $(OBJECTS) $(ZOBJECTS)
 
distclean: clean
	-rm -f $(TARGET) $(TABLEUTIL)

.PHONY: all install uninstall clean distclean
