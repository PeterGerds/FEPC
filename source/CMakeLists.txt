set(ALL_SOURCE_FILES
        basic.c
        faltung_hilfe.c
        folge.c
        fepc_cp.c
        folgen_vektor.c
        fepc_easy.c
        funktion.c
        fepc_easy_helper.c
        interval.c
        discont.c
        fepc_easy_sparse.c
        koeffizienten.c
        seconds.c
        faltung.c
        kuerzen.c)

if (HAS_FFTW3)
    set(ALL_SOURCE_FILES fft_faltung.c ${ALL_SOURCE_FILES})
endif(HAS_FFTW3)

ADD_LIBRARY(fepc SHARED ${ALL_SOURCE_FILES})

add_executable (main main.c)
target_link_libraries(main fepc m)

if (HAS_FFTW3)
    target_link_libraries(main fftw3)
endif(HAS_FFTW3)

