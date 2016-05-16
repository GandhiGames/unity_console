// from: http://virtualmore.org/wiki/index.php?title=GAContinuous
// Description
//
// Example of a continuous genetic algorithm (GA) in C# that has been ported 
// from Matlab code in "Practical Genetic Algorithms" by Haupt & Haupt, 
// 2nd Ed, (2004). A description of GAs can be found at Wikipedia.
// Usage
//
// Attach the script to a GameObject (e.g., main camera). You may then set a 
// number of parameters in the Inspector View:
//  ff --- Choose a cost function F1-F8 or F10 found in Haupt & Haupt (2004)
//  npar --- Number of optimisation variables (e.g., npar=2 for x,y)
//  varhi --- Upper limit on optimisation variables
//  varlo --- Lower limit on optimisation variables
//  maxit --- Max number of iterations
//  mincost --- Minimum cost (stop when reached)
//  popsize --- Population size
//  mutrate ---	Mutation rate
//  selection --- Fraction of population kept
//
// Press Play to observe the algorithm in action. Compare the found values with 
// the optimal ones given in Haupt & Haupt (2004).
//
// For practical use, you'll want to provide your own cost function.


using UnityEngine;
using System.Collections;
using System.Collections.Generic;

public class GAContinuous : MonoBehaviour {
	/* 
	Continuous genetic algorithm ported from Matlab version given in 
	"Practical Genetic Algorithms", Haupt & Haupt, 2nd Ed., Wiley, 2003. 
	*/
 
	// I. Setup the GA
	public int ff = 3;					// which cost function to use
	public int npar = 2; 				// number of opt vars
	public float varhi = -10; 			// Upper limit on opt vars
	public float varlo = 10;				// Lower limit on opt vars
 
	// II. Stopping criteria
	public int maxit = 100;				// max number of iterations
	public float mincost = 0; 			// minimum cost
 
	// III. GA parameters
	public int popsize = 20;		   	// set population size
	public float mutrate = 0.2f;		// set mutation rate
	public float selection = 0.5f;	   	// fraction of population kept
 
	// Display some dependent variables
	int Nt;    							// continuous parameter GA Nt=#variables
	int keep; 							// #population members that survive
	int nmut;							// total number of mutations
	int M;								// number of matings
 
	/* 
	TO DO: Make GA a coroutine with an eternal while(true) loop and do
	void Start() {
		StartCoroutine(GA());	
	}
	*/
 
	void Start() {
		GA();	
	}
 
	/* Genetic Algorithm */
	void GA() {
		/*
		Notes:
		Use PrintArray and PrintArray2D methods to examine arrays during debugging
		*/
 
		// Dependent variables
		Nt = npar;              						// continuous parameter GA Nt=#variables
		keep = Mathf.FloorToInt(selection*popsize); 	// #population members that survive
		nmut = Mathf.CeilToInt((popsize-1)*Nt*mutrate);	// total number of mutations
		M = Mathf.CeilToInt((popsize-keep)/2.0f);		// number of matings
 
		/* Create random initial population */
		float[,] par = new float[popsize,npar];
		for (int i=0; i<popsize; i++) {
			for (int j=0; j<npar; j++) {
				par[i,j] = Random.value*(varhi-varlo) + varlo;
			}
		}
 
		// Calculates initial population costs using cost function ff
		float[] costUnsorted = CostFunction(ff,par); 		
 
		// Sort the costs and obtain indices 
		float[] costSorted = new float[costUnsorted.Length];
		System.Array.Copy(costUnsorted,costSorted,costUnsorted.Length);
		System.Array.Sort(costSorted);
		int[] ind = FindIndices(costSorted, costUnsorted);	
 
		// Sort the population according to cost 
		par = SortPopulation(ind,par);
 
		// Store best (minimum) cost in a list
		List<float> minc = new List<float>();
		minc.Add(costSorted[0]);
 
		// Store average cost in a list
		List<float> meanc = new List<float>();
		meanc.Add(MeanArray(costSorted));
 
		/* Iterate through generations */
		int iga = 0;	            						// generation counter initialized
		while (iga < maxit) {
			iga++; 											// increments generation counter
 
			// Weights chromosomes
			float[] prob = new float[keep];
			float probSum = SumOfInts(1,keep);
			for (int i=0; i<keep; i++) {
				prob[i] = (i+1)/probSum;
			}
			System.Array.Reverse(prob);
 
			// Probability distribution function
			float[] odds = new float[prob.Length+1];
			float[] cumSumProb = CumSum(prob);
			for (int i=0; i<cumSumProb.Length; i++) {
				odds[i+1] = cumSumProb[i];	
			}
 
			// Mate #1 and mate #2
			float[] pick1 = new float[M];
			float[] pick2 = new float[M];
			for (int i=0; i<M; i++) {
				pick1[i] = Random.value;
				pick2[i] = Random.value;	
			}
 
			// ma and pa contain the indices of the chromosomes that will mate
			int[] ma = new int[M];
			int[] pa = new int[M];
			int ic = 0;
			while (ic<M) {
  				for (int id=1; id<keep+1; id++) {
    				if (pick1[ic]<=odds[id] && pick1[ic]>odds[id-1]) {
      					ma[ic]=id-1;
    				}
    				if (pick2[ic]<=odds[id] && pick2[ic]>odds[id-1]) {
      					pa[ic]=id-1;
    				}
  				}
  				ic++;
			}
 
			/* Performs mating using single point crossover */
			// Index of mate #1
			int keepEven = Mathf.CeilToInt(keep/2.0f);
			int[] ix = new int[keepEven]; 						
			for (int i=0; i<keepEven; i++) {
				ix[i] = i*2;	
			} 
 
			// Crossover point
			int[] xp = new int[M];
			for (int i=0; i<M; i++) {
				float rand = 0.00001f + Random.value*0.9999f; // Avoid exactly 1.0 showing up
				xp[i] = Mathf.FloorToInt(rand*Nt);	
			}
 
			// Mixing parameter
			float[] r = new float[M];
			for (int i=0; i<M; i++) {
				r[i] = Random.value;	
			}
 
			// Mating
			for (ic=0; ic<M; ic++) {				
				float xy = par[ma[ic],xp[ic]] - par[pa[ic],xp[ic]];			// ma and pa mate
				for (int i=0; i<npar; i++) {
					par[keep-1+ix[ic],i] = par[ma[ic],i];          			// 1st offspring
					par[keep+ix[ic],i] = par[pa[ic],i];          			// 2nd offspring
				}
				par[keep-1+ix[ic],xp[ic]] = par[ma[ic],xp[ic]] - r[ic]*xy; 	// 1st 
				par[keep+ix[ic],xp[ic]] = par[pa[ic],xp[ic]] + r[ic]*xy; 	// 2nd
 
				// crossover when last variable not selected
//				if (xp[ic]<npar-1) { // original if statement in Haupt not necessary?!
					for (int j=xp[ic]; j<npar-1; j++) {
						int row1 = keep-1+ix[ic];
						int row2 = keep+ix[ic];
						float tmpParRow1 = par[row1,j]; 
						par[row1,j] = par[row2,j];
//						par[row2,j] = par[row1,j]; // redundant due to line above?!	
						par[row2,j] = tmpParRow1; 
					}
//			  	}
			}
 
			/* Mutate the population */
			int[] mrow = new int[nmut];
			int[] mcol = new int[nmut];
			for (int i=0; i<nmut; i++) {
				float rand1 = 0.00001f + Random.value*0.9999f; // Avoid exactly 1.0 showing up 
				float rand2 = 0.00001f + Random.value*0.9999f; // Avoid exactly 1.0 showing up
				mrow[i] = Mathf.FloorToInt(rand1*(popsize-1))+1;	// keep row 0 (elitist strategy)
				mcol[i] = Mathf.FloorToInt(rand2*Nt);
			}
			System.Array.Sort(mrow);
 
			// Mutation
			for (int ii=0; ii<nmut; ii++) {
//			    print("mrow: " + mrow[ii]);
//				print("mcol: " + mcol[ii]);
			    par[mrow[ii],mcol[ii]] = (varhi-varlo)*Random.value+varlo;	
			}
 
			/* The new offspring and mutated chromosomes are evaluated */
			costUnsorted = CostFunction(ff,par); 						// calculate pop costs using ff
 
			/* Sort the costs and associated parameter */
			System.Array.Copy(costUnsorted,costSorted,costUnsorted.Length);
			System.Array.Sort(costSorted);
			ind = FindIndices(costSorted, costUnsorted); 
			par = SortPopulation(ind,par);
 
			/* Do statistics for a single non-averaging run */
			minc.Add(costSorted[0]);
			meanc.Add(MeanArray(costSorted));
 
			/* Stopping criteria */
//			if (iga>maxit || costSorted[0]<mincost) { // Haupt's statement is redundant!
   			if (costSorted[0]<mincost) {
   				break;
			}
 
			// Display run# and best cost so far
			print("Run: " + iga + "\t Best cost: " + costSorted[0]);
 
		} // END WHILE
 
		/* Display results */
		print("Optimized function is ff = " + ff);
		print("Best cost = " + costSorted[0]);
		PrintArrayRow("Best solution: par",par,0);	
 
	} // END START 
 
	// Cumulative sum array
	float[] CumSum(float[] array) {
		float cumSum = 0.0f;
		float[] cumArray = new float[array.Length];
		for (int i=0; i<array.Length; i++) {
			cumSum += array[i];
			cumArray[i] = cumSum;	
		}
		return cumArray;	
	}
 
	// Sum of integers from start to finish
	int SumOfInts(int start, int finish) {
		int sum = 0;
		for (int i=start; i<finish+1; i++) {
			sum += i;	
		}
		return sum;
	}
 
	// Find mean value of float list
	float MeanArray(float[] array) {
		float sumOfList = 0.0f;
		foreach (float elem in array) {
			sumOfList += elem;
		}
		return sumOfList/array.Length;
	}
 
	// Print float array with variable name aName
	void PrintArray(string aName, float[] array) {
		for (int i=0; i<array.Length; i++) {
			print(aName + "[" + i + "] = " + array[i]);
		}	
	}
 
	// Print int array with variable name aName
	void PrintArray(string aName, int[] array) {
		for (int i=0; i<array.Length; i++) {
			print(aName + "[" + i + "] = " + array[i]);
		}	
	}
 
	// Print 2D float array row by row
	void PrintArray2D(string aName, float[,] array2d) {
		int rows = array2d.GetLength(0);
		int cols = array2d.GetLength(1);
		print("cols = " + cols);
		print("rows = " + rows);
		for (int i=0; i<rows; i++) {
			string row = "";
			for (int j=0; j<cols; j++) {
				row += array2d[i,j] + " ";
			}
			print(aName + "[" + i + "] = " + row);
		}	
	}
 
	// Print row of 2D float array 
	void PrintArrayRow(string aName, float[,] array2d, int rowNumber) {
		int cols = array2d.GetLength(1);
		string row = "";
		for (int j=0; j<cols; j++) {
			row += array2d[rowNumber,j] + " ";
		}
		print(aName + "[" + rowNumber + "] = " + row);	
	}
 
 
	// Sort population in order of indices array ind
	float[,] SortPopulation(int[] ind, float[,] par) {
		int rows = par.GetLength(0);
		int cols = par.GetLength(1);
		float[,] sortedPop = new float[rows,cols];
		for (int i=0; i<rows; i++) {
			int index = ind[i];
			for (int j=0; j<cols; j++) {
				sortedPop[i,j] = par[index,j];
			}
		}		
		return sortedPop;
	} 	
 
	/* FindIndices()
	Find indices of elements in arraySorted in arrayUnsorted. Both arrays must contain
	the same set of elements. A copy of arrayUnsorted is made in tmpArray. The method
	IndexOf find the index of the *first* occurrence of an element in an array. Therefore
	that element is set to -9999 to avoid the problem of duplicates, which means that -9999
	must never be an element in the arrays. 
	*/
	int[] FindIndices(float[] arraySorted, float[] arrayUnsorted) {
		float[] tmpArray = new float[arrayUnsorted.Length];
		System.Array.Copy(arrayUnsorted,tmpArray,arrayUnsorted.Length);
		int[] indices = new int[tmpArray.Length];
		for (int i=0; i<tmpArray.Length; i++) {
			float elem = arraySorted[i];
			int index = System.Array.IndexOf(tmpArray,elem,0);
			tmpArray[index] = -9999.0f; // Requires elements to never be -9999!
			indices[i] = index;	
		}
		return indices;	
	}
 
	// Cost function to optimise over population
	float[] CostFunction(int ff, float[,] par) {
		int rows = par.GetLength(0);		
		float[] netCost = new float[rows];
 
		// TO DO: Add some cost functions 
		switch (ff) {
			// Cost function F1 (Haupt)
			case 1:
				netCost = CostF1(par);
				break;
 
			case 2:
				netCost = CostF2(par);
				break;
 
			case 3:
				netCost = CostF3(par);
				break;
 
			case 4:
				netCost = CostF4(par);
				break;
 
			case 5:
				netCost = CostF5(par);
				break;
 
			case 6:
				netCost = CostF6(par);
				break;
 
			case 7:
				netCost = CostF7(par);
				break;
 
			case 8:
				netCost = CostF8(par);
				break;
 
			case 10:
				netCost = CostF10(par);
				break;
 
			default:
				print("Cost function ff must be integer in range 1-8 or 10");
				break;
		}
 
		// Return calculated costs in an array
		return netCost;
	}
 
	// Haupt cost function F1 (p.206)
	float[] CostF1(float[,] par) {
		/*
		Requires npar=1, try varlo=-10, varhi=10
		Minimum: f(0) = 1
		*/
		int rows = par.GetLength(0);
		int cols = par.GetLength(1);		
		float[] cost = new float[rows];
		float x = new float();
		if (cols == 1) {
			for (int i=0; i<rows; i++) {
				x = par[i,0];
				cost[i] = Mathf.Abs(x) + Mathf.Cos(x);
			}
		}
		else {
			print("F1 requires npar = 1!");
		}
		return cost;
	}
 
	// Haupt cost function F2 (p.206)
	float[] CostF2(float[,] par) {
		/*
		Requires npar=1, try varlo=-100, varhi=100
		Minimum: f(0) = 0
		*/
		int rows = par.GetLength(0);
		int cols = par.GetLength(1);		
		float[] cost = new float[rows];
		float x = new float();
		if (cols == 1) {
			for (int i=0; i<rows; i++) {
				x = par[i,0];
				cost[i] = Mathf.Abs(x) + Mathf.Sin(x);
			}
		}
		else {
			print("F2 requires npar = 1!");
		}
		return cost;
	}
 
 
	// Haupt cost function F3 (p.206)
	float[] CostF3(float[,] par) {
		/*
		Try npar=2, varlo=-10, varhi=10
		Minimum: f(0,0) = 0 (note wrong value in Haupt)
 		*/
		int rows = par.GetLength(0);
		int cols = par.GetLength(1);		
		float[] cost = new float[rows];
		float x = 0.0f;
		float sum = 0.0f;
		for (int i=0; i<rows; i++) {
			sum = 0.0f;
			for (int j=0; j<cols; j++) {
				x = par[i,j];
				sum += x*x;
			}
			cost[i] = sum;
		}
		return cost;
	}
 
	// Haupt cost function F4 (p.206)
	float[] CostF4(float[,] par) {
		/*
		Try npar=2, varlo=-10, varhi=10
		Minimum: f(1,1) = 0
		*/
		int rows = par.GetLength(0);
		int cols = par.GetLength(1);		
		float[] cost = new float[rows];
		float xn = 0.0f;
		float xn1 = 0.0f;
		float sum = 0.0f;
		for (int i=0; i<rows; i++) {
			sum = 0.0f;
			for (int j=0; j<cols-1; j++) {
				xn = par[i,j];
				xn1 = par[i,j+1]; 
				sum += 100*(xn1-xn*xn)*(xn1-xn*xn)  + (1-xn)*(1-xn);
			}
			cost[i] = sum;
		}
		return cost;
	}
 
	// Haupt cost function F5 (p.207)
	float[] CostF5(float[,] par) {
		/*
		Try npar=1, varlo=-10, varhi=10
		Minimum: f(0) = -10 (note wrong value in Haupt)
		*/
		int rows = par.GetLength(0);
		int cols = par.GetLength(1);		
		float[] cost = new float[rows];
		float xn = 0.0f;
		float sum = 0.0f;
		for (int i=0; i<rows; i++) {
			sum = 0.0f;
			for (int j=0; j<cols; j++) {
				xn = par[i,j];
				sum += Mathf.Abs(xn) - 10f*Mathf.Cos(Mathf.Sqrt(Mathf.Abs(10f*xn)));
			}
			cost[i] = sum;
		}
		return cost;
	}
 
	// Haupt cost function F6 (p.207)
	float[] CostF6(float[,] par) {
		/*
		Requires npar=1, try varlo=-10, varhi=10
		Minimum: f(9.6204) = -100.22
		*/
		int rows = par.GetLength(0);
		int cols = par.GetLength(1);		
		float[] cost = new float[rows];
		float x = new float();
		if (cols == 1) {
			for (int i=0; i<rows; i++) {
				x = par[i,0];
				cost[i] = (x*x+x)*Mathf.Cos(x);
			}
		}
		else {
			print("F6 requires npar = 1!");
		}
		return cost;
	}
 
	// Haupt cost function F7 (p.207)
	float[] CostF7(float[,] par) {
		/*
		Requires npar=2, varlo=0, varhi=10
		Minimum: f(9.039,8.668) = -18.5547 (note wrong x,y in Haupt)
		*/
		int rows = par.GetLength(0);
		int cols = par.GetLength(1);		
		float[] cost = new float[rows];
		float x = 0.0f;
		float y = 0.0f;
		if (cols == 2) {
			for (int i=0; i<rows; i++) {
				x = par[i,0];
				y = par[i,1];
				cost[i] = x*Mathf.Sin(4f*x) + 1.1f*y*Mathf.Sin(2f*y);
			}
		}
		else {
			print("F7 requires npar = 2!");
		}
		return cost;
	}
 
	// Haupt cost function F8 (p.207)
	float[] CostF8(float[,] par) {
		/*
		Requires npar=2, varlo=0, varhi=10
		Minimum: f(9.039,8.668) = -18.5547 (note wrong x,y in Haupt)
		*/
		int rows = par.GetLength(0);
		int cols = par.GetLength(1);		
		float[] cost = new float[rows];
		float x = 0.0f;
		float y = 0.0f;
		if (cols == 2) {
			for (int i=0; i<rows; i++) {
				x = par[i,0];
				y = par[i,1];
				cost[i] = y*Mathf.Sin(4f*x) + 1.1f*x*Mathf.Sin(2f*y);
			}
		}
		else {
			print("F8 requires npar = 2!");
		}
		return cost;
	}
 
	// Haupt cost function F10 (p.208)
	float[] CostF10(float[,] par) {
		/*
		Try npar=2, varlo=-10, varhi=10
		Minimum: f(0,0) = 0 
		*/
		int rows = par.GetLength(0);
		int cols = par.GetLength(1);		
		float[] cost = new float[rows];
		float xn = 0.0f;
		float sum = 0.0f;
		for (int i=0; i<rows; i++) {
			sum = 0.0f;
			for (int j=0; j<cols; j++) {
				xn = par[i,j];
				sum += xn*xn - 10f*Mathf.Cos(2f*Mathf.PI*xn);
			}
			cost[i] = 10*cols + sum;
		}
		return cost;
	}
 
} // END CLASS