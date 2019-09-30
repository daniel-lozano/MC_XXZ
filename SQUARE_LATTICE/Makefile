## Pass size of the 2D array with ARGS
.PHONY:plotea

plotea:Energy_magnetization.txt plot.py
	python plot.py ${ARGS}

Energy_magnetization.txt: MC_H.out
	./MC_H.out ${ARGS}
MC_H.out:MC_XXZ.cpp
	g++ MC_XXZ.cpp -o MC_H.out

clean:
	rm Energy*.txt *.out
