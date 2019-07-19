## Pass size of the 2D array with ARGS
.PHONY:plotea

plotea:Energy_magnetization.txt plot.py
	python plot.py ${ARGS}

Energy_magnetization.txt: a.out
	./a.out ${ARGS}

a.out:MC_XXZ.cpp
	g++ MC_XXZ.cpp

clean:
	rm *.txt *.out
