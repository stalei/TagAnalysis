EXECS=ExportGalaxies
CC=gcc
all:${EXECS}

ExportGalaxies:ExportGalaxies.c
	${CC} -o ExportGalaxies ExportGalaxies.c
	
Clean:
	rm -f ${EXECS}
