CC=gcc
CFLAGS=-DHAS_FFTW3 -I/home/peet/Desktop/KO/include
LFLAGS=-L/home/peet/Desktop/KO/lib -lfftw3 -lm



all: fepcd



fepcd: main.o faltung.o kuerzen.o faltung_hilfe.o koeffizienten.o funktion.o folgen_vektor.o folge.o fft_faltung.o basic.o seconds.o
	${CC} ${CFLAGS} main.o faltung.o kuerzen.o faltung_hilfe.o koeffizienten.o funktion.o folgen_vektor.o folge.o fft_faltung.o basic.o seconds.o ${LFLAGS} -o fepcd



main.o: main.c faltung.h kuerzen.h faltung_hilfe.h koeffizienten.h funktion.h folgen_vektor.h folge.h fft_faltung.h basic.h config.h
	${CC} -c ${CFLAGS} main.c

faltung.o: faltung.c faltung.h kuerzen.h faltung_hilfe.h koeffizienten.h funktion.h folgen_vektor.h folge.h fft_faltung.h basic.h config.h
	${CC} -c ${CFLAGS} faltung.c

kuerzen.o: kuerzen.c kuerzen.h faltung_hilfe.h koeffizienten.h funktion.h folgen_vektor.h folge.h fft_faltung.h basic.h config.h
	${CC} -c ${CFLAGS} kuerzen.c

faltung_hilfe.o: faltung_hilfe.c faltung_hilfe.h koeffizienten.h funktion.h folgen_vektor.h folge.h fft_faltung.h basic.h config.h
	${CC} -c ${CFLAGS} faltung_hilfe.c

koeffizienten.o: koeffizienten.c koeffizienten.h funktion.h folgen_vektor.h folge.h fft_faltung.h basic.h config.h
	${CC} -c ${CFLAGS} koeffizienten.c

funktion.o: funktion.c funktion.h folgen_vektor.h folge.h fft_faltung.h basic.h config.h
	${CC} -c ${CFLAGS} funktion.c

folgen_vektor.o: folgen_vektor.c folgen_vektor.h folge.h fft_faltung.h basic.h config.h
	${CC} -c ${CFLAGS} folgen_vektor.c

folge.o: folge.c folge.h fft_faltung.h basic.h config.h
	${CC} -c ${CFLAGS} folge.c

fft_faltung.o: fft_faltung.c fft_faltung.h basic.h config.h
	${CC} -c ${CFLAGS} fft_faltung.c

basic.o: basic.c basic.h config.h
	${CC} -c ${CFLAGS} basic.c

seconds.o: seconds.c seconds.h
	${CC} -c ${CFLAGS} seconds.c



clean:
	rm -f *~ *.o fepcd



