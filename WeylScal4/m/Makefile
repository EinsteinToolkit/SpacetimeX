# -*-Makefile-*-

all: WeylScal4.out
	@echo
	@echo "The Cactus thorns are up to date."
	@echo

WeylScal4.out: WeylScal4.m
	rm -rf WeylScal4
	./runmath.sh $^
	grep -q 'Thorn ./WeylScal4 created successfully' WeylScal4.out && \
	for file in $$(ls WeylScal4); do				\
		./copy-if-changed.sh WeylScal4/$$file ../$$file;	\
	done

clean:
	rm -rf WeylScal4
	rm -f WeylScal4.out WeylScal4.err

.PHONY: all clean
