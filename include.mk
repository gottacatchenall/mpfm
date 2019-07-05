ROOT_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

DIR_PREFIX:=$(ROOT_DIR)/src/core

MAIN =              $(DIR_PREFIX)/main.cpp

INCLUDE_DIRS = 		$(DIR_PREFIX)						\
					$(DIR_PREFIX)/MPFM_Instance					\
					$(DIR_PREFIX)/Population			\
					$(DIR_PREFIX)/DataWrangler			\
					$(DIR_PREFIX)/Individual

SRCS =				$(DIR_PREFIX)/MPFM_Instance/MPFM_Instance.cpp					\
					$(DIR_PREFIX)/Population/Population.cpp		\
					$(DIR_PREFIX)/DataWrangler/DataWrangler.cpp 	\
					$(DIR_PREFIX)/Individual/Individual.cpp
