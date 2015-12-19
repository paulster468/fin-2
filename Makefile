#
# Download and place this file into your git repository. 
# Rename it to 'Makefile' (without quotes)
#
# wget http://www.phys.uconn.edu/phys2200/downloads/Makefile.sample
# mv Makefile.sample Makefile
#
# or
#
# wget http://www.phys.uconn.edu/phys2200/downloads/Makefile.sample -O Makefile
#

EDITOR = /usr/bin/vim
INDENT = /usr/bin/indent

CC        = clang
LDFLAGS   = -O
CFLAGS    = -Wall -Wextra -pedantic $(LDFLAGS)
LDLIBS    = $(shell gsl-config --libs)
GNUPLOT = gnuplot

.SUFFIXES:
.SUFFIXES:  .c .o .h

.PHONY: edit clean veryclean

target = mainf

$(target) : $(target).o two-cubes-integrand2.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

edit : $(target).c
	$(EDITOR) $<
	$(INDENT) $<

res: mainf
	./mainf > res

show: res cubes.gp
	$(GNUPLOT) cubes.gp
	eog cubes.png

clean : 
	rm -f *.o
	rm -f *.*~

veryclean : clean
	rm  -f $(target)

