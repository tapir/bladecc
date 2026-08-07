CC      ?= gcc
CFLAGS  ?= -O2 -Wall -Wextra -std=c99
CPPFLAGS += -Iinclude
LDLIBS  := -lm

LIB     := libbladecc.a
LIBSRC  := src/interp.c src/airfoils.c src/ccblade.c
LIBOBJ  := $(LIBSRC:.c=.o)

CPPCHECK       ?= cppcheck
CPPCHECK_FLAGS ?= --enable=all --inconclusive --std=c99 --force --quiet \
                  -Iinclude --inline-suppr \
                  --suppress=missingIncludeSystem \
                  --suppress=checkersReport \
                  --suppress=normalCheckLevelMaxBranches \
                  --error-exitcode=1

SCAN_BUILD      ?= scan-build
SCAN_BUILD_FLAGS ?= --status-bugs

CLANG_FORMAT      ?= clang-format
FORMAT_SOURCES    := include/ccblade.h src/interp.c src/airfoils.c \
                     src/ccblade.c tests/runtests.c tests/compare_c.c

.PHONY: all test clean check analyze format format-check

all: $(LIB)

$(LIB): $(LIBOBJ)
	ar rcs $@ $^

%.o: %.c include/ccblade.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

tests/runtests: tests/runtests.c $(LIB) include/ccblade.h tests/reference_data.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ tests/runtests.c $(LIB) $(LDLIBS)

test: tests/runtests
	cd tests && ./runtests

check:
	$(CPPCHECK) $(CPPCHECK_FLAGS) include/ src/ tests/runtests.c tests/compare_c.c

analyze: clean
	$(SCAN_BUILD) $(SCAN_BUILD_FLAGS) $(MAKE) tests/runtests

format:
	$(CLANG_FORMAT) -i $(FORMAT_SOURCES)

format-check:
	$(CLANG_FORMAT) --dry-run --Werror $(FORMAT_SOURCES)

clean:
	rm -f $(LIBOBJ) $(LIB) tests/runtests
