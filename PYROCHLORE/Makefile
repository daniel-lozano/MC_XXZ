.PHONY:plotea

plotea: Energy_${ARGx}_${ARGy}_${ARGz}_H.txt  plot_closed_loop.py
	python plot_closed_loop.py ${ARGx} ${ARGy} ${ARGz}

Energy_${ARGx}_${ARGy}_${ARGz}_H.txt : loop.out
	./loop.out 

loop.out:MC_3NN.cpp 
	g++ MC_3NN.cpp -o loop.out

clean:
	rm *.out
	rm *.txt
