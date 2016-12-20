
import java.util.*;
import java.io.*;
//import org.apache.commons.math3.*;
//args syntax: simulated #ofsamples mutationList.txt numclusters trials

public class ConsensusTreeBuilding {
    public static class Tree {
        ArrayList<ArrayList<String>> clusters;
        boolean[][] adjacencyMatrix;
        double[] mutationalFrequency;
        int numClusters;
        double[][] distanceMatrix;
        int coverage = 100;

		public Tree(ArrayList<ArrayList<String>> clusters, boolean[][] adjacencyMatrix, int numClusters, double[] frequencies, double[][] distanceMatrix) {
			this.clusters = clusters;
            this.adjacencyMatrix = adjacencyMatrix;
            this.numClusters = numClusters;
            mutationalFrequency = frequencies;
            this.distanceMatrix = distanceMatrix;
		}

		public ArrayList<ArrayList<String>> getClusters() {
			return clusters;
		}

		public boolean[][] getAdjacencyMatrix() {
			return adjacencyMatrix;
		}
        
        public int getNumClusters() {
            return numClusters;
        }
        
        public double[] getFrequencies() {
            return mutationalFrequency;
        }
        
        public double[][] getDistanceMatrix() {
            return distanceMatrix;
        }
        
        public int getCoverage() {
            return coverage;
        }
        
        public void setCoverage(int c) {
            coverage = c;
        }
	}
     
    public static void main(String[] args) { 
        ArrayList<String> mutationList = new ArrayList<String>();
        ArrayList<Tree> testData = new ArrayList<Tree>();
        int numSamples = Integer.parseInt(args[1]);
        int trials = 1; //for tests on real data
        int numClusters = 5; //if no number of clusters entered
        if(args.length==5) {
            trials = Integer.parseInt(args[4]);
            numClusters =Integer.parseInt(args[3]);;
        }
        double[][] printMatrix = new double[trials][104];
        for(int k=0; k<trials; k++) {
            System.out.println(k);
            ArrayList<Integer> pr = new ArrayList<Integer>();
            //read in the trees from a dataset
            if(args[0].equals("real")) {
                mutationList = commonMutations(args[2]);
                ArrayList<Tree> trees = readAllDataFile(args[2], mutationList);
                //trees before they are reduced to only include mutations that
                //are included in all trees
                for(int i=0; i<trees.size(); i++) {
                    String treeName = "originalInput" + i +".dot";
                    writeDOT(trees.get(i).getAdjacencyMatrix(), treeName, trees.get(i).getClusters(), pr);
                }
                testData = reduceTrees(trees, mutationList);
            } 


            if(args[0].equals("simulated")) {
                String mutationFile = args[2];
                mutationList = readMutation(mutationFile);
                mutationList.add(0, "Normal Cell");
            }

            //randomly create 'true' tree
            ArrayList<Integer> differenceList = new ArrayList<Integer>();
            ArrayList<Double> freq = new ArrayList<Double>();
            ArrayList<ArrayList<String>> clusters = createClusters(numClusters, mutationList);
            boolean [][] adjacencyMatrix = makeMatrix(numClusters, clusters);
            double [] frequencies = makeFrequencies(adjacencyMatrix, clusters, mutationList);
            double[][] distanceMatrix = makeDistanceMatrix(mutationList, clusters, adjacencyMatrix);
            Tree trueTree = new Tree(clusters, adjacencyMatrix, numClusters, frequencies, distanceMatrix);
            ArrayList<Tree> trueList = new ArrayList<Tree>();    
            trueList.add(trueTree);

            //create test data set of size numSamples, with some error
            if(args[0].equals("simulated")) {
                testData = createTestData(trueTree, numSamples, mutationList);
                //draw true tree
                writeDOT(adjacencyMatrix, "trueTreefromTest.dot", clusters, pr);
            }

            //draw the input trees
            for(int i=0; i<numSamples; i++) {
                String treeName = "sampleinput" + i +".dot";
                writeDOT(testData.get(i).getAdjacencyMatrix(), treeName, testData.get(i).getClusters(), pr);
            }

            //creating majority vote tree

            //create matrix of incidences of parent of each mutation
            int[][] complexParentMatrix = createParentMatrix(testData, mutationList);
            int[][] complexClusterMatrix = createClusterMatrix(testData, mutationList);
            double majority = numSamples/(2.0); //threshold
            int[][] parentMatrix = parseMatrix(complexParentMatrix, majority);
            int[][] clusterMatrix = parseMatrix(complexClusterMatrix, majority);
            int[][] childMatrix = createChildMatrix(parentMatrix);

            //build the new clusters and adjacency matrix
            ArrayList<ArrayList<Integer>> newClusters = buildNewClusters(parentMatrix, childMatrix, clusterMatrix, complexClusterMatrix);
            boolean[][] newAdjacencyMatrix = buildNewMatrix(newClusters, parentMatrix);

            //tree before conflicts are resolved:
            ArrayList<ArrayList<String>> newClusterNames = createClustersNames(newClusters, mutationList);
            int[] frequencyCheckMatrix = frequencyCheck(testData, newClusterNames, newAdjacencyMatrix, mutationList);
            ArrayList<Integer> problems = printProblemClusters(frequencyCheckMatrix, newClusterNames, mutationList, numSamples);
            //writeDOT(newAdjacencyMatrix, "MajorityVoteConflict.dot", newClusterNames, problems); 
            frequencyCheckMatrix = frequencyCheck(testData, newClusterNames, newAdjacencyMatrix, mutationList);

            //Majority vote tree after conflicts are resolved:
            ArrayList<Integer> conflicts = findConflicts(newClusters, newAdjacencyMatrix, mutationList);
            newAdjacencyMatrix = resolveConflicts(newClusters, newAdjacencyMatrix, conflicts, frequencyCheckMatrix, mutationList);
            newAdjacencyMatrix = connectTree(newAdjacencyMatrix, newClusters, complexParentMatrix, mutationList);
            newClusterNames = createClustersNames(newClusters, mutationList);
            double[] notReal = makeFrequencies(newAdjacencyMatrix, newClusterNames, mutationList);
            double[][] newDistanceMatrix = makeDistanceMatrix(mutationList, newClusterNames, newAdjacencyMatrix);
            Tree majorityVoteTree =  new Tree(newClusterNames, newAdjacencyMatrix, newClusterNames.size(), notReal, newDistanceMatrix);
            ArrayList<Tree> finalTree = new ArrayList<Tree>();
            finalTree.add(majorityVoteTree);
            int[][] consensusMatrix = createParentMatrix(finalTree, mutationList);
            problems = printProblemClusters(frequencyCheckMatrix, newClusterNames, mutationList, numSamples);
            writeDOT(newAdjacencyMatrix, "MajorityVoteTreeNow.dot", newClusterNames, pr); 


           //tree created w/ Distance-Based Markov Chain
            Tree mcDistanceTree = runMCMCDistance(testData, mutationList, majorityVoteTree);
            mcDistanceTree = reroot(mcDistanceTree);
            frequencyCheckMatrix = frequencyCheck(testData, mcDistanceTree.getClusters(), mcDistanceTree.getAdjacencyMatrix(), mutationList);
            boolean[][] a1 = ancestor(mcDistanceTree.getAdjacencyMatrix(), mutationList, mcDistanceTree.getClusters());

            problems = printProblemClusters(frequencyCheckMatrix, mcDistanceTree.getClusters(), mutationList, numSamples);
            writeDOT(mcDistanceTree.getAdjacencyMatrix(), "mcTreeDistanceTreeNow.dot", mcDistanceTree.getClusters(), pr);

            //tree created with New Metric and MCMC process
            double alpha = .2;
            Tree mcmcTree1 = runMCMC(testData, mutationList, majorityVoteTree, 0);
            mcmcTree1 = reroot(mcmcTree1);
            writeDOT(mcmcTree1.getAdjacencyMatrix(), "mcmcTree1.dot", mcmcTree1.getClusters(), pr);

            Tree mcmcTree2 = runMCMC(testData, mutationList, majorityVoteTree, .3);
            mcmcTree2 = reroot(mcmcTree2);
            writeDOT(mcmcTree2.getAdjacencyMatrix(), "mcmcTree2.dot", mcmcTree2.getClusters(), pr);

            Tree mcmcTree3 = runMCMC(testData, mutationList, majorityVoteTree, .5);
            mcmcTree3 = reroot(mcmcTree3);
            writeDOT(mcmcTree3.getAdjacencyMatrix(), "mcmcTree3.dot", mcmcTree3.getClusters(), pr);

            Tree mcmcTree4 = runMCMC(testData, mutationList, majorityVoteTree, .7);
            mcmcTree4 = reroot(mcmcTree4);
            writeDOT(mcmcTree4.getAdjacencyMatrix(), "mcmcTree4.dot", mcmcTree4.getClusters(), pr);

            Tree mcmcTree5 = runMCMC(testData, mutationList, majorityVoteTree, 1);
            mcmcTree5 = reroot(mcmcTree5);
            writeDOT(mcmcTree5.getAdjacencyMatrix(), "mcmcTree5.dot", mcmcTree5.getClusters(), pr);

            // boolean[][] a2 = ancestor(mcmcTree.getAdjacencyMatrix(), mutationList, mcmcTree.getClusters());
            // boolean[][] a = ancestor(trueTree.getAdjacencyMatrix(), mutationList, trueTree.getClusters());

            //comparison between true tree and consensus trees and input trees
            if(args[0].equals("simulated")) {
                printMatrix = outputMatrix(trueTree, majorityVoteTree, mcDistanceTree, testData, mutationList, k, printMatrix, mcmcTree1, mcmcTree2, mcmcTree3, mcmcTree4, mcmcTree5);
            }
            
        
        }
        //write comparisons/evaluations to csv file: outputFile.csv
        if(args[0].equals("simulated")) {
                printEvaluations(printMatrix);
            }

    }
    
    //reads in the list of mutations and adds them to an ArrayList, which is returned
    public static ArrayList<String> readMutation(String mutationFile) {
        ArrayList<String> mutationList = new ArrayList<String>();
        try {
			File f = new File(mutationFile);
			Scanner sc = new Scanner(f);
            String curLine = sc.nextLine();
            while (sc.hasNextLine()) {
				curLine = sc.nextLine();
                mutationList.add(curLine);
            }
        }
            catch(Exception ex) {
				ex.printStackTrace();
		}
        return mutationList;
    }
    
    
    //randomly adds each mutation to a cluster
    public static ArrayList<ArrayList<String>> createClusters(int numClusters, ArrayList<String> mutationList) {
        ArrayList<ArrayList<String>> clusters = new ArrayList<ArrayList<String>>();
        for (int i=0; i<numClusters+1; i++) {
            clusters.add(new ArrayList<String>());
        }
        clusters.get(0).add(mutationList.get(0));
        Random rn = new Random();
        //adds a mutation to each cluster
        for(int i=1; i<numClusters+1; i++) {
            clusters.get(i).add(mutationList.get(i));
        }
        //add the rests of the mutations to a random cluster
        for (int i=numClusters+1; i<mutationList.size(); i++) {
            int location =rn.nextInt(numClusters)+1;
            clusters.get(location).add(mutationList.get(i));
        }
        return clusters;    
    }
    
    //sets cluster 0 as the root, the iterates through the remaining clusters 
    //and randomly select an existing node as the parent for that cluster
    //returns an adjacency matrix of the tree
    public static boolean[][] makeMatrix(int numClusters, ArrayList<ArrayList<String>> clusters) {
        boolean[][] adjacencyMatrix = new boolean[numClusters+1][numClusters+1];
        ArrayList<Integer> addedNodes = new ArrayList<Integer>();
        int root = 0;
        addedNodes.add(root);
        adjacencyMatrix[0][1] = true;
        addedNodes.add(1);
        for(int y=2; y<numClusters+1; y++) {
            int x = (int)(Math.random() * addedNodes.size()-1)+1;
            adjacencyMatrix[x][y] = true;
            addedNodes.add(y);
        }
        return adjacencyMatrix;
    }
    
    //creates the frequencies for a tree - sum of the frequencies of children must be
    //less than the frequency of the parent
    public static double[] makeFrequencies(boolean[][] adjacencyMatrix, ArrayList<ArrayList<String>> clusters, ArrayList<String> mutationList) {
        double[] frequencies = new double[mutationList.size()];
        ArrayList<Integer> parents = new ArrayList<Integer>();
        int root = 0;
        parents.add(root);
        frequencies[root] = 100.0;
        int i =0;
        while(i < adjacencyMatrix.length && i<parents.size()) {
            double upperBound = frequencies[parents.get(i)];
            for(int j =0; j<adjacencyMatrix.length; j++) {
                if(adjacencyMatrix[parents.get(i)][j]) {
                    double lowerBound = upperBound/5.0;
                    if(lowerBound ==0) {
                        lowerBound = .01;
                    }
                    double childF = (Math.random()*(upperBound-lowerBound)) + lowerBound;
                    if((upperBound-childF)<1) {
                        upperBound =1;
                    } else {
                        upperBound = upperBound - childF;
                    }
                    for(String mutation:clusters.get(j)) {
                        int index = mutationList.indexOf(mutation);
                        frequencies[index] = childF;
                    }
                    parents.add(j);
                }
            }
            i++;
        }
        return frequencies;   
        
    }
    
    //returns a random sample from a Binomial Distribution with parameters n and p
    public static int getBinomial(int n, double p) {
        int x = 0;
        Random r = new Random();
        for(int i = 0; i < n; i++) {
            double rand = r.nextDouble();
            if(rand < p)
                x++;
        } 
        return x;
        }

    //returns a random sample from a Poisson Distribution centered at lambda
    public static int getPoisson(double lambda) {
        double L = Math.exp(-lambda);
        double p = 1.0;
        int k = 0;
        Random r = new Random();
        while (p > L) {
            k++;
            p *= r.nextDouble();
        }
        return k - 1;
    }
    
     //creates dataset of trees with some error
    public static ArrayList<Tree> createTestData(Tree trueTree, int numSamples, ArrayList<String> mutationList) {
        ArrayList<Tree> testData = new ArrayList<Tree>();
        for(int i=0; i<numSamples; i++) {
            Tree sample = createTestTree(trueTree, mutationList);
            testData.add(sample);
        }
        return testData;
    }

    //creates a tree with some error
    public static Tree createTestTree(Tree trueTree, ArrayList<String> mutationList){
        boolean[][] tempAdjacencyMatrix = trueTree.getAdjacencyMatrix();
        ArrayList<ArrayList<String>> testClusters = new ArrayList<ArrayList<String>>();
        int numClusters = trueTree.getNumClusters();
        boolean[][] testAdjacencyMatrix = new boolean[numClusters+1][numClusters+1];
        double[] frequencies = new double[mutationList.size()];
        double[] trueFrequencies = trueTree.getFrequencies();
        //create adjacency matrix for new tree -> starts as the same as the true tree
        for(int x=0; x<numClusters+1; x++) {
            for(int y=0; y<numClusters+1; y++) {
                testAdjacencyMatrix[x][y] = tempAdjacencyMatrix[x][y];
            }
        }
        //create clusters -> starts as the same as the true tree
        for(int x=0; x<numClusters+1; x++) {
            ArrayList<String> testGroup = new ArrayList<String>();
            for(int y=0; y<trueTree.getClusters().get(x).size(); y++) {
                testGroup.add(trueTree.getClusters().get(x).get(y));
            }
            testClusters.add(testGroup);
        }
        
        //create new frequencies (I use numbers between 0 and 100, but the actual mutational frequencies are that number /100)
        int coverage = getPoisson(100);
        frequencies[0] = 100.0;
        for(int x=1; x<frequencies.length; x++) {
            double adjustedF = trueFrequencies[x];
            double newFrequency = trueFrequencies[x];
            Random r = new Random();
            double rand = r.nextDouble();
            //only change some proportion of the frequencies
            if(rand<.3) {
                double p = adjustedF/100.0;
                int variantAlleleReads = getBinomial(coverage, p);
                double referenceReads = (double) coverage - variantAlleleReads;
                double f = (variantAlleleReads)/(variantAlleleReads + referenceReads);
                newFrequency = f*100;
                if(newFrequency ==0.0) {
                    newFrequency =.001;
                }
                if(newFrequency >= 100) {
                    newFrequency = 99;
                }
            }
            frequencies[x] = newFrequency;
        }

        
        //find errors that correspond to changes in frequencies and change tree accordingly - step 1
        ArrayList<String> outliers = new ArrayList<String>();
        ArrayList<Integer> clusterChanges = new ArrayList<Integer>();
        ArrayList<Double> averageFrequencies = calculateAverageFrequencies(testClusters, frequencies, mutationList);
        
        for(int x=1; x<tempAdjacencyMatrix.length; x++) {
            for(String m:testClusters.get(x)) {
                int ind = mutationList.indexOf(m);
                //if a mutation has a mutational frequency that is 10% away from the average of that cluster (the population frequency),
                //then it should not be in the cluser
                double lowerBound = Math.min((averageFrequencies.get(x)*.9), averageFrequencies.get(x)-1);
                double upperBound = Math.max((averageFrequencies.get(x)*1.1), averageFrequencies.get(x)+1);
                if(((frequencies[ind]<lowerBound) || frequencies[ind]>upperBound) && testClusters.get(x).size()>1){   
                    outliers.add(m);
                    clusterChanges.add(x);
                    }
            }
        }

        for(int i=0; i<outliers.size(); i++) {
            //see if its frequency is similar to another cluster, if so move it to that cluster
            int newCluster = -1;
            for(int j=1; j<averageFrequencies.size(); j++) {
                if(frequencies[mutationList.indexOf(outliers.get(i))] > averageFrequencies.get(j)*0.9 && frequencies[mutationList.indexOf(outliers.get(i))] < averageFrequencies.get(j)*1.1) {
                    newCluster = j;
                }
            }
            if(newCluster != -1) {
                //move outliers.get(i) to newCluster - error 3
                if(clusterChanges.get(i) != newCluster) {
                    
                    if(testClusters.get(clusterChanges.get(i)).size() >1) {
                        
                        testClusters = moveMutations(testClusters, testClusters.size(), testAdjacencyMatrix, clusterChanges.get(i), newCluster, outliers.get(i));
                        averageFrequencies = calculateAverageFrequencies(testClusters, frequencies, mutationList);
                    }
                    
                }
             } else {
                //if not, move it to be parent or child
                if(averageFrequencies.get(clusterChanges.get(i)) > frequencies[mutationList.indexOf(outliers.get(i))]) {
                    //move outliers.get(i) to be a child - error 4
                    if(testClusters.get(clusterChanges.get(i)).size()>1) {
                    
                        testAdjacencyMatrix = addChild(testClusters, testAdjacencyMatrix, clusterChanges.get(i), outliers.get(i));

                        averageFrequencies = calculateAverageFrequencies(testClusters, frequencies, mutationList);
                    }
                } else {
                    //move outliers.get(i) to be parent - error 4
                    if(testClusters.get(clusterChanges.get(i)).size()>1) {
                        testAdjacencyMatrix = addParent(testClusters, testAdjacencyMatrix, clusterChanges.get(i), outliers.get(i));
                       averageFrequencies = calculateAverageFrequencies(testClusters, frequencies, mutationList);
                    }
                }
            }
        }
        
        //step 2 -> check that parent frequency is greater than children frequency
        ArrayList<Integer> next = new ArrayList<Integer>();
        next.add(0);
        int current =0;
        while(current < testAdjacencyMatrix.length && current<next.size()) {
            for(int i=0; i<testAdjacencyMatrix.length; i++) {
                if(testAdjacencyMatrix[next.get(current)][i]) {
                    if(averageFrequencies.get(i) > averageFrequencies.get(next.get(current)) && (next.get(current) !=0)) {
                        if(averageFrequencies.get(next.get(current))<(averageFrequencies.get(i)*.8)) {
                                //move i to be the new parent -> error 2 (.8 is arbitrary cutoff)
                                testAdjacencyMatrix = moveNodes(testAdjacencyMatrix, testAdjacencyMatrix.length, (next.get(current)), i);
                                averageFrequencies = calculateAverageFrequencies(testClusters, frequencies, mutationList);
                            } else{
                                //switch the two clusters -> error1
                                testAdjacencyMatrix = mutateNodes(tempAdjacencyMatrix, testAdjacencyMatrix, numClusters, (next.get(current)), i);
                                averageFrequencies = calculateAverageFrequencies(testClusters, frequencies, mutationList);
                            }
                        next.add(next.get(current));
                        next.set(current, i);
                    } else {
                        next.add(i);
                    }
                }
            }
            current++;
        }
        
        //step 3 --> check for sum rule 
        ArrayList<Integer> n = new ArrayList<Integer>();
        n.add(0);
        int counter = 0;
        while(counter < testAdjacencyMatrix.length && counter<n.size()) {
            boolean again = true;
            while(again) { 
                again=false;
                ArrayList<Integer> children = new ArrayList<Integer>();
                double childSum = 0;
                double lowest = 100.0;
                int lowestIndex = 0;
                double highest = 0.0;
                int highestIndex = 0;
                for(int i=0; i<testAdjacencyMatrix.length; i++) {
                    if(testAdjacencyMatrix[n.get(counter)][i]) {
                        childSum += averageFrequencies.get(i);
                       if(averageFrequencies.get(i) > highest) {
                            highest = averageFrequencies.get(i);
                            highestIndex = i;
                        }
                        if(averageFrequencies.get(i) < lowest) {
                            lowest = averageFrequencies.get(i);
                            lowestIndex = i;
                        }
                        n.add(i);
                        children.add(i);
                    }
                }
                if(childSum > averageFrequencies.get(n.get(counter))) {
                    //move child - move clusters with the lowest frequency to be child of cluster with the highest frequency
                    //repeatedly move clusters until it adheres to the sum rule
                    double difference = childSum - averageFrequencies.get(n.get(counter));
                    if(difference>0 && lowest<difference) {
                        again=true;
                    }
                    if(lowestIndex != highestIndex) {
                        testAdjacencyMatrix = moveNodes(testAdjacencyMatrix, testAdjacencyMatrix.length, lowestIndex, highestIndex);
                        averageFrequencies = calculateAverageFrequencies(testClusters, frequencies, mutationList);
                    }
                    //make sure highest comes before lowest in next (need list to be ordered so that a 
                    //clusters parent is always earlier in the list)
                    if(n.indexOf((Integer)lowestIndex) < n.indexOf((Integer)highestIndex)) {
                        int temp = n.indexOf((Integer)highestIndex);
                        int temp2 = n.indexOf((Integer)lowestIndex);
                        n.set(temp, lowestIndex);
                        n.set(temp2, highestIndex);
                    }

                }
            }
            counter++;
        }
        //step 4 --> combine a similar cluster
        boolean combine = true;
        for(int i=1; i<averageFrequencies.size()-1; i++) {
            for(int j=i+1; j<averageFrequencies.size(); j++) {
                if(combine && testAdjacencyMatrix[i][j] && (averageFrequencies.get(i) <= averageFrequencies.get(j)*1.1)) {
                    //combine clusters i and j
                    testAdjacencyMatrix = combineClusters(testClusters, testAdjacencyMatrix, i, j);
                    averageFrequencies = calculateAverageFrequencies(testClusters, frequencies, mutationList);
                    combine = false;
                } if(combine && testAdjacencyMatrix[j][i] && averageFrequencies.get(i)*1.1 >= averageFrequencies.get(j)) {
                    //combine clusters i and j
                    testAdjacencyMatrix = combineClusters(testClusters, testAdjacencyMatrix, j, i);
                    averageFrequencies = calculateAverageFrequencies(testClusters, frequencies, mutationList);
                    combine = false;
                }
            }
        }

        double[][] distanceMatrix = makeDistanceMatrix(mutationList, testClusters, testAdjacencyMatrix);   
        Tree testTree = new Tree(testClusters, testAdjacencyMatrix, numClusters, frequencies, distanceMatrix);
        testTree.setCoverage(coverage);
        return testTree;
        

    }
    
    //calculates the population frequencies for each cluster, which is equal to the average of the mutational frequencies
    //of all mutations in that cluster
    public static ArrayList<Double> calculateAverageFrequencies(ArrayList<ArrayList<String>> testClusters, double[] frequencies, ArrayList<String> mutationList) {
        ArrayList<Double> averageFrequencies = new ArrayList<Double>();
        averageFrequencies.add(100.0);
        for(int x=1; x<testClusters.size(); x++) {
            double averageFrequency = 0.0;
            for(String mutation: testClusters.get(x)) {
                int index = mutationList.indexOf(mutation);
                averageFrequency += frequencies[index];
            }
            averageFrequency = averageFrequency/testClusters.get(x).size();
            averageFrequencies.add(averageFrequency);
        }
        return averageFrequencies;
    }

    //error1: switch parent and child node in the tree
    public static boolean[][] mutateNodes(boolean[][] tempAdjacencyMatrix, boolean[][] testAdjacencyMatrix, int numClusters, int cluster1, int cluster2) {
        for (int i=0; i<testAdjacencyMatrix.length; i++) {
            boolean tempValue1 = testAdjacencyMatrix[cluster1][i];
            testAdjacencyMatrix[cluster1][i] = testAdjacencyMatrix[cluster2][i];
            testAdjacencyMatrix[cluster2][i] = tempValue1;
        }
          
        for (int i=0; i<testAdjacencyMatrix.length; i++) {
            boolean tempValue = testAdjacencyMatrix[i][cluster1];
            testAdjacencyMatrix[i][cluster1] = testAdjacencyMatrix[i][cluster2];
            testAdjacencyMatrix[i][cluster2] = tempValue;
        }
        return testAdjacencyMatrix;
    }
    
    //error 2: make a child (or sibling) node the parent of its parent (or sibling) in the tree
    public static boolean[][] moveNodes(boolean[][] testAdjacencyMatrix, int numClusters, int node1, int newParent) {
        int oldParent = node1;
        for(int i=0; i<testAdjacencyMatrix.length; i++) {
            if(testAdjacencyMatrix[i][newParent]) {
                oldParent = i;
            }
        }
        //parent is not longer parent of child
        testAdjacencyMatrix[oldParent][newParent] = false;
        
        //parent of parent is now parent of child    
        for (int i=0; i<numClusters; i++) {
            if(testAdjacencyMatrix[i][node1]) {
                testAdjacencyMatrix[i][newParent] = true;
                testAdjacencyMatrix[i][node1] = false;
            }
        }        
        //child becomes parent of parent
        testAdjacencyMatrix[newParent][node1] = true;
        
        return testAdjacencyMatrix;
    }
    
    //error3: move mutations from one cluster to the parent or child of the cluster
    public static ArrayList<ArrayList<String>> moveMutations(ArrayList<ArrayList<String>> testClusters, int numClusters, boolean[][] testAdjacencyMatrix, int cluster1, int cluster2, String m) {
        Random r = new Random();
        if(testClusters.get(cluster1).size()>1) {
            String moveMutation = m;
            int mutation = 0;
            if(m.equals("")) {
                mutation = r.nextInt(testClusters.get(cluster1).size());
                moveMutation = testClusters.get(cluster1).get(mutation);   
            } else {
                mutation = testClusters.get(cluster1).indexOf(moveMutation);
            }
            testClusters.get(cluster1).remove(mutation);
            testClusters.get(cluster2).add(moveMutation);
    
        }
        return testClusters;
    }
    //error 4ish: take a mutation from one cluster and create a new cluster at a random location
    public static boolean[][] addCluster(ArrayList<ArrayList<String>> testClusters, boolean[][] testAdjacencyMatrix, int node1, String outlier) {
        boolean [][] largerAdjacencyMatrix = new boolean[testAdjacencyMatrix.length + 1][testAdjacencyMatrix.length + 1];
        for(int i =0; i<testAdjacencyMatrix.length; i++) {
            for(int j=0; j<testAdjacencyMatrix.length; j++) {
                largerAdjacencyMatrix[i][j] = testAdjacencyMatrix[i][j];
            }
        }
        int clusterLoc = node1;
        int mutation = testClusters.get(clusterLoc).indexOf(outlier);
        String newMutation = testClusters.get(clusterLoc).get(mutation);
        int newCluster = (int)(Math.random() * testAdjacencyMatrix.length);
        largerAdjacencyMatrix[newCluster][largerAdjacencyMatrix.length-1] = true;
        testClusters.get(clusterLoc).remove(mutation);
        ArrayList<String> newAddedCluster = new ArrayList<String>();
        newAddedCluster.add(newMutation);
        testClusters.add(newAddedCluster);
        
        return largerAdjacencyMatrix;
    }
    
    //error 5: combine a child and parent cluster
    public static boolean[][] combineClusters(ArrayList<ArrayList<String>> testClusters, boolean[][] testAdjacencyMatrix, int node1, int node2) {
        boolean [][] smallerAdjacencyMatrix = new boolean[testAdjacencyMatrix.length -1][testAdjacencyMatrix.length - 1];
        int parent = node1;
        int child = node2;
        if(parent!=child) {
        //add mutations from child to parent
        for(int i=0; i<testClusters.get(child).size(); i++) {
            testClusters.get(parent).add(testClusters.get(child).get(i));
        }
        testClusters.remove(child);
        //add children of child to parent
        for (int i=0; i<testAdjacencyMatrix.length; i++) {
            if (testAdjacencyMatrix[child][i]) {
                testAdjacencyMatrix[parent][i] = true;
            }
        }
        int xIndex =0;
        //re fill in table
        for (int i=0; i<smallerAdjacencyMatrix.length; i++) {
            if(xIndex == child) {
                    xIndex++;
                }
            int yIndex =0;
            for (int j=0; j<smallerAdjacencyMatrix.length; j++) {
                if(yIndex == child) {
                    yIndex++;
                }
                smallerAdjacencyMatrix[i][j] = testAdjacencyMatrix[xIndex][yIndex];
                yIndex++;
            }
            xIndex++;
        }
        }
        return smallerAdjacencyMatrix;
    }
    
   
    
    //creates a matrix that marks the number of times mutation x is the parent of mutation y
    public static int[][] createParentMatrix(ArrayList<Tree> testTrees, ArrayList<String> mutationList) {
        int [][] parentMatrix = new int[mutationList.size()][mutationList.size()];
        //for each test tree
        for(int i=0; i<testTrees.size(); i++) {
            ArrayList<ArrayList<String>> sampleClusters = testTrees.get(i).getClusters();
            boolean[][] matrix = testTrees.get(i).getAdjacencyMatrix();
            //x=the cluster
            for(int x=0; x<matrix.length; x++) {
                //y=each mutation in the cluster
                for(int y=0; y<sampleClusters.get(x).size(); y++) {
                    //index of mutation in parentMatrix
                    int index = mutationList.indexOf(sampleClusters.get(x).get(y));
                    //look at all clusters
                    for(int z=0; z<matrix.length; z++) {
                        //if cluster x is the parent of cluster z
                        if(matrix[x][z]) {
                            //increment every mutation in cluster z in the matrix by one
                            for(int w=0; w<sampleClusters.get(z).size(); w++) {
                                int childIndex = mutationList.indexOf(sampleClusters.get(z).get(w));
                                parentMatrix[index][childIndex] = (parentMatrix[index][childIndex]) + 1;
                               
                            }
                        }
                    }
                }
            }
        }
        return parentMatrix;
    }
    
    //creates a matrix that marks the number of times mutation y is in the same cluster as mutation z
    public static int[][] createClusterMatrix(ArrayList<Tree> testTrees, ArrayList<String> mutationList) {
        int [][] clusterMatrix = new int[mutationList.size()][mutationList.size()];
        //for each test tree
        for(int i=0; i<testTrees.size(); i++) {
            ArrayList<ArrayList<String>> sampleClusters = testTrees.get(i).getClusters();
            //x=the cluster
            for(int x=0; x<sampleClusters.size(); x++) {
                //y=each mutation in the cluster
                for(int y=0; y<sampleClusters.get(x).size(); y++) {
                    //index of mutation y
                    int index = mutationList.indexOf(sampleClusters.get(x).get(y));
                    //z=other mutations in the cluster
                    for(int z=0; z<sampleClusters.get(x).size(); z++) {
                        int otherIndex = mutationList.indexOf(sampleClusters.get(x).get(z));
                        //if mutation y and z are not the same
                        if(y!=z) {
                            clusterMatrix[index][otherIndex] = (clusterMatrix[index][otherIndex]) + 1;
                               
                            }
                        }
                    }
                }
            }
        
        return clusterMatrix;
    }
    

    //if matrix[x][y] > threshold, then mutation x is the parent or sibling (depending on input matrix)
    //of mutation y in the majority of trees
    //If it meets the threshold, then the value is a 1, otherwise it is a 0
    public static int[][] parseMatrix(int[][] matrix, double majority) {
        int[][] simplifiedMatrix = new int[matrix.length][matrix.length];
        for(int i=0; i<matrix.length; i++) {
            for(int j=0; j<matrix[0].length; j++) {
                if(matrix[i][j] > majority) {
                    simplifiedMatrix[i][j] = 1;
                }
                else {
                    simplifiedMatrix[i][j] = 0;
                }
            }  
        } 
        return simplifiedMatrix;
    }
    
    //uses parentMatrix to create childMatrix
    public static int[][] createChildMatrix(int[][] parentMatrix) {
        int[][] childMatrix = new int[parentMatrix.length][parentMatrix[0].length];
        for(int i=0; i<parentMatrix.length; i++) {
            for(int j=0; j<parentMatrix[0].length; j++) {
                childMatrix[i][j] = parentMatrix[j][i];
            }
        }
        return childMatrix;
    }
    
    //determine clustering of mutations for the majority vote tree
    public static ArrayList<ArrayList<Integer>> buildNewClusters(int[][] parentMatrix, int[][] childMatrix, int[][] clusterMatrix, int[][] complexClusterMatrix) {
        ArrayList<ArrayList<Integer>> newClusters = new ArrayList<ArrayList<Integer>>();
        ArrayList<Integer> alreadyPaired = new ArrayList<Integer>();
        for(int i=0; i<parentMatrix.length; i++) {
            ArrayList<Integer> pairings = new ArrayList<Integer>();
            if(!alreadyPaired.contains(i)){
                pairings.add(i);
                alreadyPaired.add(i);
                for(int j=0; j<parentMatrix[0].length; j++) {
                    //if two mutations have the same parent, children and are siblings in the majority of trees, then they are clustered together
                    if (Arrays.equals(parentMatrix[i], parentMatrix[j]) && Arrays.equals(childMatrix[i], childMatrix[j]) && i!=j && clusterMatrix[i][j] ==1) {
                        pairings.add(j);
                        alreadyPaired.add(j);
                    }
                }
            newClusters.add(pairings);
            }
        }
        for(int x=0; x<parentMatrix.length; x++) {
            if(!alreadyPaired.contains(x)) {
                int mostSimilar = 0;
                int similarityScore = 0;
                for(int i=0; i<clusterMatrix.length; i++) {
                    if(complexClusterMatrix[x][i]>similarityScore) {
                        similarityScore = complexClusterMatrix[x][i];
                        mostSimilar = i;
                    }
                }
                for(int j=0; j<newClusters.size(); j++) {
                    if(newClusters.get(j).contains(mostSimilar)) {
                        newClusters.get(j).add(x);
                        alreadyPaired.add(x);
                    }
                }
                if(!alreadyPaired.contains(x)) {
                    ArrayList<Integer> addCluster = new ArrayList<Integer>();
                    addCluster.add(x);
                    newClusters.add(addCluster);
                    alreadyPaired.add(x);
                }
                
            }
        }
        return newClusters;
    }
    
    //use new clusterings to create new adjacency matrix
    public static boolean[][] buildNewMatrix(ArrayList<ArrayList<Integer>> newClusters, int[][] parentMatrix) {
        boolean[][] newMatrix = new boolean[newClusters.size()][newClusters.size()];
        for(int i=0; i<newClusters.size(); i++) {
            for(int j=0; j<newClusters.size(); j++) {
                for(int x=0; x<newClusters.get(i).size(); x++) {
                    for(int y=0; y<newClusters.get(j).size(); y++) {
                        //any edge between two mutations that occurs in the majority of trees, is added to the new tree
                        if(newClusters.get(i).size()>0 && newClusters.get(j).size()>0 && parentMatrix[newClusters.get(i).get(x)][newClusters.get(j).get(y)] == 1) {
                            newMatrix[i][j] = true;
                        }
                    }
                }
            }
        }
        return newMatrix;
    }
    
    //find clusters that are in conflict with the restrictions of a tumor phylogeny
    //ex. two cluster share the same children, two clusters contain the same mutation, cycles
    public static ArrayList<Integer> findConflicts(ArrayList<ArrayList<Integer>> newClusters, boolean[][] newMatrix, ArrayList<String> mutationList) {
        ArrayList<Integer> inConflict = new ArrayList<Integer>();
        ArrayList<Integer> inCluster = new ArrayList<Integer>();
        ArrayList<ArrayList<String>> clusterNames = createClustersNames(newClusters, mutationList);
        boolean [][] ancestors = ancestor(newMatrix, mutationList, clusterNames);

        
        for(int i=0; i<newMatrix.length; i++) {
              if(newMatrix[i][0]) {
                newMatrix[i][0] = false;
              }
        }
        for(int i=0; i<newMatrix.length; i++) {
            for(int j=0; j<newMatrix.length; j++) {
                for(int k=0; k<newMatrix.length; k++) {
                    //if cluster i and cluster j share a child or point at each other/are part of a cycle, then they are in conflict
                    int index1 = newClusters.get(i).get(0);
                    int index2 = newClusters.get(j).get(0);
                    if(((newMatrix[i][k] && newMatrix[j][k]) || (newMatrix[i][j] && newMatrix[j][i]) || (ancestors[index1][index2] && ancestors[index2][index1])) && i!=j) {
                        if(!inConflict.contains(i)) {
                            inConflict.add(i);
                        } 
                        if(!inConflict.contains(j)) {
                            inConflict.add(j);

                        }
                    }
                }
                
            }    
        }
        //Adds to conflict if two clusters share the same mutation 
        for(int i=0; i<newClusters.size(); i++) {
            for(int j=0; j<newClusters.get(i).size(); j++) {
                if(inCluster.contains(newClusters.get(i).get(j))) {
                   if(!inConflict.contains(i)) {
                        inConflict.add(i);
                   }
                   for(int c=0; c<newClusters.size(); c++) {
                        if(newClusters.get(c).contains(newClusters.get(i).get(j)) && c!=i && !inConflict.contains(c)) {
                            inConflict.add(c);   
                        }
                   }
                
                }
                inCluster.add(newClusters.get(i).get(j));
            }
        }
        return inConflict;
    }
    
    //merge clusters that share the same children
    public static boolean[][] resolveConflicts(ArrayList<ArrayList<Integer>> newClusters, boolean[][] newMatrix, ArrayList<Integer> inConflict, int[] frequencyCheckMatrix, ArrayList<String> mutationList) {
        ArrayList<Integer> combinedCluster = new ArrayList<Integer>();
        ArrayList<ArrayList<String>> clusterNames = createClustersNames(newClusters, mutationList);
        boolean [][] ancestors = ancestor(newMatrix, mutationList, clusterNames);
        if(inConflict.size()>0) {
            int mostCommon = -1; 
            int numShared = 0;
            int numDifferent = 0;
            int sameParent = 0;
            int conflict1 = inConflict.get(0);
            ArrayList<Integer> conflictChildren = newClusters.get(conflict1);
            //find the cluster that conflict1 is most similar to 
            //(aka shares the most children with, is part of the same cycle or share a mutation)
            for(Integer c:inConflict) {
                int currentDifferent = 0;
                int currentShared = 0;
                int parentShared = 0;
                for(int j=0; j<newMatrix.length; j++) {
                    if(newMatrix[conflict1][j] && newMatrix[c][j]) {
                        currentShared += 1;
                    } else if(newMatrix[conflict1][j] || newMatrix[c][j]) {
                        currentDifferent +=1;
                    } 
                }
                int index1 = newClusters.get(conflict1).get(0);
                int index2 = newClusters.get(c).get(0);
                //cycle
                if((newMatrix[conflict1][c] && newMatrix[c][conflict1]) || (ancestors[index1][index2] && ancestors[index2][index1])) {
                    currentShared +=100;
                    currentDifferent = 0;
                }
                
    
                //repeat mutation
                for(int m: newClusters.get(conflict1)) {
                        if(newClusters.get(c).contains(m) && c!=conflict1) {
                            currentShared +=100;
                            currentDifferent = 0;
                        }
                }
                
                if(currentShared > numShared && c!=conflict1) {
                        numShared = currentShared;
                        mostCommon = c;
                        numDifferent = currentDifferent;
                        sameParent = parentShared;
                       
                }
                
                
            }
            //combine mutations in the similar clusters to create a single cluster
            if(numDifferent<=2 && conflict1!=0 && mostCommon!=0 && numShared !=0) {
                
                boolean [][] smallerAdjacencyMatrix = new boolean[newMatrix.length -1][newMatrix.length - 1];
                int temp = mostCommon;
                if(mostCommon < conflict1) {
                    mostCommon = conflict1;
                    conflict1 = temp;
                }
                for(int k =0; k<newClusters.get(mostCommon).size(); k++) {
                    combinedCluster.add(newClusters.get(mostCommon).get(k));
                }
                for(int k =0; k<newClusters.get(conflict1).size(); k++) {
                    if(!combinedCluster.contains(newClusters.get(conflict1).get(k))){
                        combinedCluster.add(newClusters.get(conflict1).get(k));
                    }
                }
                newClusters.set(conflict1, combinedCluster);
                newClusters.remove(mostCommon);

                //add children and parents of cluster2 to cluster1, and get rid of cluster2                                     
                int xIndex =0;
                for(int x=0; x<smallerAdjacencyMatrix.length; x++) {
                    if(xIndex == mostCommon) {
                        xIndex++;
                    }
                    int yIndex =0;
                    for (int j=0; j<smallerAdjacencyMatrix.length; j++) {
                        if(yIndex == mostCommon) {
                            yIndex++;
                        }
                        smallerAdjacencyMatrix[x][j] = newMatrix[xIndex][yIndex];
                        if(newMatrix[mostCommon][yIndex]) {
                            smallerAdjacencyMatrix[conflict1][j] = true;
                        }
                        yIndex++;
                        }

                    if(newMatrix[xIndex][mostCommon]) { 
                        smallerAdjacencyMatrix[x][conflict1] = true;
                    }
                    xIndex++;  
                }
                //make sure cluster1 isn't a parent to itself
                smallerAdjacencyMatrix[conflict1][conflict1] = false;
                newMatrix = smallerAdjacencyMatrix;

               //check for remaining conflicts
                inConflict = findConflicts(newClusters, newMatrix, mutationList);

            }
            //remove an edge (if clusters are not similar, but share a child)
            else {
                int child = -1;
                int remove = conflict1;
                int frequencySum1 = 0;
                int frequencySum2 = 0;
                for(int m:newClusters.get(conflict1)) {
                    frequencySum1 += frequencyCheckMatrix[m];
                }
                for(int m:newClusters.get(mostCommon)) {
                    frequencySum2 += frequencyCheckMatrix[m];
                }
                frequencySum1 = frequencySum1/newClusters.get(conflict1).size();
                frequencySum2 = frequencySum2/newClusters.get(mostCommon).size();
                if(frequencySum2 < frequencySum1) {
                    remove = mostCommon;
                }
                for(int x=0; x<newMatrix.length; x++) {
                    if(newMatrix[conflict1][x] && newMatrix[mostCommon][x]){
                        child = x;
                    }
                }
                
                
             
                newMatrix[remove][child] = false;
                
                //check for remaining conflicts
                inConflict = findConflicts(newClusters, newMatrix, mutationList);
            }
        } else {
            
            return newMatrix;
        }
        //repeatedly resolve conflicts until none remain
        return resolveConflicts(newClusters, newMatrix, inConflict, frequencyCheckMatrix, mutationList);

    }
    
    
    //add most likely parent to any cluster that is not the normal cell and does not have a parent
    //turns forest into a tree
    public static boolean[][] connectTree(boolean[][] newMatrix, ArrayList<ArrayList<Integer>> newClusters,  int[][] parentMatrix, ArrayList<String> mutationList) {
        ArrayList<ArrayList<String>> clusterNames = createClustersNames(newClusters, mutationList);
        boolean [][] ancestors = ancestor(newMatrix, mutationList, clusterNames);
        for(int i=1; i<newMatrix.length; i++) {
            boolean noParent = true;
            for(int j=0; j<newMatrix.length; j++) {
                if(newMatrix[j][i]) {
                    noParent = false;
                }
            }
            int mostLikely = 0;
            int likelihood = 0;
            int index1 = 0;
            int index2 =0;
            if(noParent) {
                for(int x=0; x<newClusters.size(); x++) {
                    int currentLikelihood = 0;
                    for(int c2: newClusters.get(x)) {
                        for(int c1: newClusters.get(i)) {
                            currentLikelihood += parentMatrix[c2][c1];
                        }
                    }
                    index1 = newClusters.get(i).get(0);
                    index2 = newClusters.get(x).get(0);
                    currentLikelihood = currentLikelihood/newClusters.get(x).size();
                    //cluster can't be the parent of cluster i if it is a descendant of i (no cycles)
                    if(currentLikelihood>likelihood && x!=i && !ancestors[index2][index1]) {
                        mostLikely = x;
                        likelihood = currentLikelihood;
                    }
                }
                newMatrix[mostLikely][i] = true;
                ancestors = ancestor(newMatrix, mutationList, clusterNames);
            }
        }
        return newMatrix;
    }
    
    //returns all the clusters which are not well supported based on the mutational frequencies
    public static ArrayList<Integer> printProblemClusters(int[] frequencyCheck, ArrayList<ArrayList<String>> newClusters, ArrayList<String> mutationList, int numSamples) {
        ArrayList<Integer> problems = new ArrayList<Integer>();
        double sum = 0.0;
        for(int i=0; i<frequencyCheck.length; i++) {
            for(int j=0; j<newClusters.size(); j++) {
                String mutation = mutationList.get(i);
                if(newClusters.get(j).contains(mutation)) {
                    double problem = frequencyCheck[i]/((double)numSamples);
                    if(problem <= .5) {
                        problems.add(j);
                    }
                    sum += problem;
                }
            }
        }
        return problems;
    }


    //change id to actual name of mutation in list of clusters
    public static ArrayList<ArrayList<String>> createClustersNames(ArrayList<ArrayList<Integer>> newClusters, ArrayList<String> mutationList){
        ArrayList<ArrayList<String>> clusterNames = new ArrayList<ArrayList<String>>();
        for(int i=0; i<newClusters.size(); i++) {
            ArrayList<String> groupNames = new ArrayList<String>();
            for(int j=0; j<newClusters.get(i).size(); j++) {
                String mutationName = mutationList.get(newClusters.get(i).get(j));
                groupNames.add(mutationName);
            }
            clusterNames.add(groupNames);
        }
        return clusterNames;
    }
    
    //makes a distance matrix for a given tree (distance between any two mutations in the tree)
    public static double[][] makeDistanceMatrix(ArrayList<String> mutationList, ArrayList<ArrayList<String>> clusters, boolean[][] adjacencyMatrix) {

        double[][] distanceMatrix = new double[mutationList.size()][mutationList.size()];
        int root=0;
        ArrayList<Integer> children = new ArrayList<Integer>();
        children.add(root);
        int index = 0;
        int numberOfClusters = clusters.size();
        while(index < numberOfClusters) {
            ArrayList<String> currentCluster = clusters.get(children.get(index));
            double d = 0;
            //distance of 0 for all mutations in the same cluster
            for(String m:currentCluster) {
                for(String m2:currentCluster) {
                    int location = mutationList.indexOf(m);
                    int location2 = mutationList.indexOf(m2);
                    distanceMatrix[location][location2] = d;
                    }
            }
            d++;
            ArrayList<Integer> next = new ArrayList<Integer>();
            next.add(children.get(index));
            //all descendents
            while(next.size() > 0) {
                ArrayList<Integer> current = new ArrayList<Integer>();
                for(int x=0; x<next.size(); x++) {
                    int currentIndex = next.get(x);
                    for(int i=0; i<adjacencyMatrix.length; i++) {
                        if(adjacencyMatrix[currentIndex][i]) {
                            if(index==0) {
                                children.add(i);
                            }
                            //children of the current one
                            current.add(i);
                            ArrayList<String> otherCluster = clusters.get(i);
                            for(String m:currentCluster) {
                                for(String m2:otherCluster) {
                                    int location = mutationList.indexOf(m);
                                    int location2 = mutationList.indexOf(m2);
                                    distanceMatrix[location][location2] = d;
                                    }
                            }
                        }
                    }
                }
                next=current;
                //increment distance
                d++;
            }
            //mutation m will have a distance to mutations that are not descendants of m, that is
            //one greater than the distance between m's parent and that mutation 
            int parent = -1;
            for(int y=0; y<adjacencyMatrix.length; y++) {
                if(adjacencyMatrix[y][children.get(index)]) {
                    parent = y;
                }
            }
            if(parent != -1 && clusters.get(parent).size() > 0) {
                int parentIndex = mutationList.indexOf(clusters.get(parent).get(0));
                for(String mLoc:currentCluster) {
                    int mIndex = mutationList.indexOf(mLoc);
                    for(int k=0; k<mutationList.size(); k++) {
                        if(distanceMatrix[mIndex][k] == 0 && !currentCluster.contains(mutationList.get(k))) {
                            distanceMatrix[mIndex][k] = distanceMatrix[parentIndex][k] +1;
                        }
                    }
                }
            }
            
            if(index==0) {
               numberOfClusters = children.size();
            }
            index++;
        }
        return distanceMatrix;
    }
    
   
    
    //returns the distance between two trees
    public static int calculateDistance(double[][] distanceMatrix1, double[][] distanceMatrix2) {
        int difference = 0;
        for(int i=0; i<distanceMatrix1.length-1; i++) {
            for(int j=i+1; j<distanceMatrix1[0].length; j++) {
                int currentDifference = (int) Math.abs(distanceMatrix1[i][j] - distanceMatrix2[i][j]);
                difference += currentDifference;
            }
        }
        return difference;
    }
    
    //checks whether the structure of the consensus tree is possible given all the frequencies seen from the test data
    public static int[] frequencyCheck(ArrayList<Tree> testTrees, ArrayList<ArrayList<String>> clusters, boolean[][] adjacencyMatrix, ArrayList<String> mutationList) {
    //for each parent child relationship in parentMatrix check if it is possible given the frequencies from each tree in test trees
    //for each tree
        int[] isPossible = new int[mutationList.size()];
        for(Tree t:testTrees) {
            double[] frequencies = t.getFrequencies();
            for(int i =0; i<mutationList.size(); i++) {
                double parentFrequency = frequencies[i];
                double childFrequency = 0.0;
                for(int j =0; j<adjacencyMatrix.length; j++) {
                    int index=0;
                    for(int x=0; x<adjacencyMatrix.length; x++) {
                        if(clusters.get(x).contains(mutationList.get(i))) {
                            index = x;
                        }
                    }
                    if(adjacencyMatrix[index][j]) {
                        double clusterFrequency = 0.0;
                        for(String mutation:clusters.get(j)) {
                            int ind = mutationList.indexOf(mutation);
                            clusterFrequency += frequencies[ind];
                        }
                        clusterFrequency = clusterFrequency/clusters.get(j).size();
                        childFrequency += clusterFrequency;
                    }  
                }
                if(parentFrequency >= childFrequency) {
                    isPossible[i] = isPossible[i] + 1;   
                } 
            }
        }
        return isPossible;
    }
    
   

   
    //calculates the total distance between a tree and all the test trees
    public static int distanceSum(ArrayList<Tree> testTrees, Tree currentTree) {
        int sum = 0;
        for(Tree t: testTrees) {
            int distance = calculateDistance(t.getDistanceMatrix(), currentTree.getDistanceMatrix());
            sum = sum + distance;
        }
        return sum;
    }
        
    //creates a 'neighbor' tree to the current tree
    //a cell in the distance matrix for a neighbor tree differs
    //by no more than 1 value from a cell in the distance matrix of the current tree
    public static Tree createNeighbor(Tree currentTree, ArrayList<String> mutationList) {
        boolean[][] currentAdjacencyMatrix = currentTree.getAdjacencyMatrix();
        ArrayList<ArrayList<String>> newClusters = new ArrayList<ArrayList<String>>();
        int numClusters = currentAdjacencyMatrix.length;
        boolean[][] newAdjacencyMatrix = new boolean[numClusters][numClusters];
        double[] frequencies = new double[mutationList.size()];
        //create adjacency matrix for new tree
        for(int x=0; x<numClusters; x++) {
            for(int y=0; y<numClusters; y++) {
                newAdjacencyMatrix[x][y] = currentAdjacencyMatrix[x][y];
            }
        }
        //create clusters
        for(int x=0; x<numClusters; x++) {
            ArrayList<String> group = new ArrayList<String>();
            for(int y=0; y<currentTree.getClusters().get(x).size(); y++) {
                group.add(currentTree.getClusters().get(x).get(y));
            }
            newClusters.add(group);
        }
        Random rand = new Random();

        
        //choose change
        int change = rand.nextInt(9);
        //find a random parent and child cluster
        int parent = 0;
        int child = -1;
        int sibling =-1;
        while(child==-1) {
            parent = rand.nextInt(newClusters.size());
            ArrayList<Integer> children = new ArrayList<Integer>();
            ArrayList<Integer> parents = new ArrayList<Integer>();
            ArrayList<Integer> siblings = new ArrayList<Integer>();
            int parentOf = 0;
            for(int j=0; j<newAdjacencyMatrix.length; j++) {
                if(newAdjacencyMatrix[j][parent]) {
                    parentOf = j;
                }
            }
            for(int i=0; i<newAdjacencyMatrix.length; i++) {
                if(newAdjacencyMatrix[parent][i]) {
                    children.add(i);
                }
                if(newAdjacencyMatrix[parentOf][i] && i!=parent) {
                    siblings.add(i);
                }
                
            }
            if(siblings.size() > 0) {
                int r = rand.nextInt(siblings.size());
                sibling = siblings.get(r);
            }
            if(children.size() > 0) {
                int r = rand.nextInt(children.size());
                child = children.get(r);
            } else {
                for(int i=0; i<newAdjacencyMatrix.length; i++) {
                    if(newAdjacencyMatrix[i][parent]) {
                        parents.add(i);
                    }
                }
                if(parents.size() >0) {
                    int rn = rand.nextInt(parents.size());
                    child = parent;
                    parent = parents.get(rn);
                }
            }
        }
        //if it has no sibling
        if(change == 8 && sibling == -1) {
            change = 1;
        }
        //if the tree only has one cluster containing mutations
        if(change==2 && newAdjacencyMatrix.length ==2) {
            change = 7;
        }
        
        
        //find a random cluster
        int randomCluster = rand.nextInt(newClusters.size());
        int randomCluster2 = rand.nextInt(newClusters.size());
        while(randomCluster == randomCluster2) {
            randomCluster2 = rand.nextInt(newClusters.size());
        }
       
        
        
        //switch parent and child
        if(change ==0) {
            newAdjacencyMatrix = mutateNodes(currentAdjacencyMatrix, newAdjacencyMatrix, newAdjacencyMatrix.length, parent, child);  
        //move child to be parent of parent
        }else if (change==1) {
            newAdjacencyMatrix = moveNodes(newAdjacencyMatrix, newAdjacencyMatrix.length, parent, child);
        //combine clusters
        }else if (change==2) {
            newAdjacencyMatrix = combineClusters(newClusters, newAdjacencyMatrix, parent, child);
        //make child the sibling of parent
        }else if (change==3) {
         newAdjacencyMatrix = makeSiblings(newAdjacencyMatrix, newAdjacencyMatrix.length, parent, child);   
        //move mutation from child to parent
        }else if (change==4) {
            if(newClusters.get(child).size()>1) {
                newClusters = moveMutations(newClusters, newAdjacencyMatrix.length, newAdjacencyMatrix, child, parent, "");
            }
        //move mutation from parent to child
        }else if (change==5) {
            if(newClusters.get(parent).size()>1) {
                newClusters = moveMutations(newClusters, newAdjacencyMatrix.length, newAdjacencyMatrix, parent, child, "");
            }
        //move mutation from random cluster to new cluster that is parent    
        }else if (change==6) {
            if(newClusters.get(randomCluster).size() >1) {
                newAdjacencyMatrix = addParent(newClusters, newAdjacencyMatrix, randomCluster, "");
                
            }
        //move mutation from random cluster to new cluster that is a child
        }else if (change==7) {
            if(newClusters.get(randomCluster).size() >1) {
                newAdjacencyMatrix = addChild(newClusters, newAdjacencyMatrix, randomCluster, "");
            }
         
        } else {
            //sibling becomes parent of its sibling
            newAdjacencyMatrix = moveNodes(newAdjacencyMatrix, newAdjacencyMatrix.length, parent, sibling);

        } 
        double[][] distanceMatrix = makeDistanceMatrix(mutationList, newClusters, newAdjacencyMatrix);
        Tree neighbor = new Tree(newClusters, newAdjacencyMatrix, newAdjacencyMatrix.length, frequencies, distanceMatrix);
        return neighbor;
                
    }
    
    //takes a mutation and makes a new cluster that is child to the original cluster
    public static boolean[][] addChild(ArrayList<ArrayList<String>> testClusters, boolean[][] testAdjacencyMatrix, int parent, String outlier) {
        boolean [][] largerAdjacencyMatrix = new boolean[testAdjacencyMatrix.length + 1][testAdjacencyMatrix.length + 1];
        if(testClusters.get(parent).size() > 0) {
        for(int i =0; i<testAdjacencyMatrix.length; i++) {
            for(int j=0; j<testAdjacencyMatrix.length; j++) {
                largerAdjacencyMatrix[i][j] = testAdjacencyMatrix[i][j];
            }
        }
        String newMutation = outlier;
        int mutation = 0;
        if(newMutation.equals("")) {   
            Random r = new Random();
            mutation = r.nextInt(testClusters.get(parent).size());
            newMutation = testClusters.get(parent).get(mutation);
        } else {
            mutation = testClusters.get(parent).indexOf(newMutation);
        }
        largerAdjacencyMatrix[parent][largerAdjacencyMatrix.length-1] = true;
        testClusters.get(parent).remove(mutation);
        ArrayList<String> newAddedCluster = new ArrayList<String>();
        newAddedCluster.add(newMutation);
        testClusters.add(newAddedCluster);
        }
        return largerAdjacencyMatrix;
    }
    
    //takes a mutation and makes a new cluster that is parent to the original cluster
    public static boolean[][] addParent(ArrayList<ArrayList<String>> testClusters, boolean[][] testAdjacencyMatrix, int cluster, String outlier) {
        boolean [][] largerAdjacencyMatrix = new boolean[testAdjacencyMatrix.length + 1][testAdjacencyMatrix.length + 1];
        for(int i =0; i<testAdjacencyMatrix.length; i++) {
            for(int j=0; j<testAdjacencyMatrix.length; j++) {
                largerAdjacencyMatrix[i][j] = testAdjacencyMatrix[i][j];
            }
        }
        String mutation = outlier;
        int random = 0;
        if(mutation.equals("")) {
            random = (int) Math.random() * testClusters.get(cluster).size();
            mutation = testClusters.get(cluster).get(random);
        } else {
            random = testClusters.get(cluster).indexOf(mutation);
        }
        int parent = cluster;
        int child = largerAdjacencyMatrix.length -1;
        //parent of parent is now parent of child    
        for (int i=0; i<testAdjacencyMatrix.length; i++) {
            if(testAdjacencyMatrix[i][parent]) {
                largerAdjacencyMatrix[i][child] = true;
                largerAdjacencyMatrix[i][parent] = false;
            }
        }        
        //child becomes parent of parent
        largerAdjacencyMatrix[child][parent] = true;
        
        ArrayList<String> addCluster = new ArrayList<String>();
        addCluster.add(mutation);
        testClusters.add(addCluster);
        testClusters.get(parent).remove(random);
        return largerAdjacencyMatrix;
    }
    
    //takes child and makes it the sibling of its parent
    public static boolean[][] makeSiblings(boolean[][] testAdjacencyMatrix, int numClusters, int parent, int child) {
        //parent of parent is now parent of child 
        int oldParent = parent;
        for(int i=0; i<testAdjacencyMatrix.length; i++) {
            if(testAdjacencyMatrix[i][child]) {
                oldParent = i;   
            }
        }
        for (int i=0; i<numClusters; i++) {
            if(testAdjacencyMatrix[i][parent]) {
                testAdjacencyMatrix[i][child] = true;
                //parent is not longer parent of child        
                testAdjacencyMatrix[oldParent][child] = false;
            }
        }        
        return testAdjacencyMatrix;
    }
    
    
    //creates an initial random tree to begin with when running the MCMC
    public static Tree createInitialTree(ArrayList<String> mutationList) {
        //currently must manually enter a starting number of clusters (although the value doesn't affect it too much)
        int numClusters = 5;
        ArrayList<ArrayList<String>> clusters = createClusters(numClusters, mutationList);
        boolean [][] adjacencyMatrix = makeMatrix(numClusters, clusters);
        double[][] distanceMatrix = makeDistanceMatrix(mutationList, clusters, adjacencyMatrix);
        double[] frequencies = new double[mutationList.size()];
        Tree initialTree = new Tree(clusters, adjacencyMatrix, numClusters, frequencies, distanceMatrix);
        return initialTree;
    }

    
    //uses a markov chain monte carlo to find the minimum complete-data distance tree
    //the complete-data distance of a tree is found by calculating the distance between that tree
    //and all other trees from the test data
    public static Tree runMCMCDistance(ArrayList<Tree> testTrees, ArrayList<String> mutationList, Tree initial) {
        ArrayList<Integer> problems = new ArrayList<Integer>();
        //can randomly initialize the tree, or initialize it with the majority vote tree
       // Tree current = createInitialTree(mutationList);
       // Tree heated = createInitialTree(mutationList);
       // Tree hot = createInitialTree(mutationList);
       // Tree hotter = createInitialTree(mutationList);

        //4 different chains of the tree
        //only uses neighboring tree
        Tree current = initial;
        //uses 2-distance neighbors 
        Tree heated = initial;
        //uses 3-distance neighbors 
        Tree hot = initial;
        //uses 4-distance neighbors 
        Tree hotter = initial;

        int p1 = 0;
        int i=0;
        double min = 9999999999.0;
        Tree best = current;
        while(i<10000) {
            Random r = new Random();
            p1 = distanceSum(testTrees, current);
            Tree neighbor = createNeighbor(current, mutationList);
            int p2 = distanceSum(testTrees, neighbor);
            double random = r.nextDouble();
            double ratio = ((double)p1)/p2;
            int heatedP = distanceSum(testTrees, heated);
            Tree heatedNeighbor = createNeighbor(heated, mutationList);
            int heatedP2 = distanceSum(testTrees, heatedNeighbor);
            double heatedRatio = ((double)heatedP)/heatedP2;
            int hotP = distanceSum(testTrees, hot);
            Tree hotNeighbor = createNeighbor(hot, mutationList);
            hotNeighbor = createNeighbor(hotNeighbor, mutationList);
            int hotP2 = distanceSum(testTrees, hotNeighbor);
            double hotRatio = ((double)hotP)/hotP2;
            int hotterP = distanceSum(testTrees, hotter);
            Tree hotterNeighbor = createNeighbor(hotter, mutationList);
            int hotterP2 = distanceSum(testTrees, hotterNeighbor);
            double hotterRatio = ((double)hotterP2)/hotterP;
            if(p2<p1) {
                current = neighbor;
                p1=p2;
            }
            double random2 = r.nextDouble();
            if(heatedP2<heatedP) {
                heated = heatedNeighbor;
                heatedP = heatedP2;
            }
            double random3 = r.nextDouble();
             if(hotP2<hotP) {
                hot = hotNeighbor;
                hotP = hotP2;
            }
            double random4 = r.nextDouble();
             if(hotterP2<hotterP) {
                hotter = hotterNeighbor;
                hotterP = hotterP2;
            }
            if(heatedP<p1) {
                Tree temp = current;
                current = heated;
                heated = temp;
                p1 = heatedP;
            }
            if(hotP<p1) {
                Tree temp = current;
                current = hot;
                hot = temp;
                p1 = hotP;
            }
            if(hotterP<p1) {
                Tree temp = current;
                current = hotter;
                hotter = temp;
                p1 = hotterP;
            }
            if(p1<min) {
                min=p1;
                best=current;
            }
            i++;
            
        }
        
        return current;
    }
    
    //uses a markov chain monte carlo to find the tree the the highest ancestral similarity and lowest distance
 //to the trees in the test data
    public static Tree runMCMC(ArrayList<Tree> testTrees, ArrayList<String> mutationList, Tree initial, double alpha) {
        ArrayList<Integer> problems = new ArrayList<Integer>();
        // Tree current = createInitialTree(mutationList);
        // Tree heated = createInitialTree(mutationList);
        Tree current = initial;
        Tree heated = initial;
        //metric: alpha(x) + (1-alpha)y, where x is the ancestral probability and y is the distance probability
        double p1 = 0.0;
        int i=0;
        double min = 9999999999.0;
        Tree best = current;
        while(i<10000) {
            Random r = new Random();
            Tree neighbor = createNeighbor(current, mutationList);
            p1 = combinedComparison(testTrees, current, mutationList, alpha);
            double p2 = combinedComparison(testTrees, neighbor, mutationList, alpha);
            double random = r.nextDouble();
            double ratio = p2/p1;
            double heatedP = combinedComparison(testTrees, heated, mutationList, alpha);
            Tree heatedNeighbor = createNeighbor(heated, mutationList);
            double heatedP2 = combinedComparison(testTrees, heatedNeighbor, mutationList, alpha);
            double heatedRatio = heatedP/heatedP2;
            double random2 = r.nextDouble();
            if(p2>p1 /*|| (random<ratio)*/) {
                current = neighbor;
                p1=p2;
            }
            if(heatedP2>heatedP /* || (random2<heatedRatio)*/) {
                heated = heatedNeighbor;
                heatedP = heatedP2;
            }
            if(heatedP>p1) {
                Tree temp = current;
                current = heated;
                p1=heatedP;
                heated = current;
            }
            if(p1>min) {
                min=p1;
                best=current;
            }
            i++;
        }
        
        return current;
    }
    
    //a combined measure to evaluate the similarity in terms of similarity in ancestral relationships
    //and distance
    public static double combinedComparison(ArrayList<Tree> testData, Tree tree, ArrayList<String> mutationList, double alpha) {
        double ancestralAverage = 0.0;
        double distanceAverage = 0.0;
        boolean[][] a = ancestor(tree.getAdjacencyMatrix(), mutationList, tree.getClusters());
        for(Tree t:testData) {
            boolean[][] currentA = ancestor(t.getAdjacencyMatrix(), mutationList, t.getClusters());
            double currentAncestral = ancestralDifferences(a, currentA);
            ancestralAverage += currentAncestral;
            int currentDistance = calculateDistance(tree.getDistanceMatrix(), t.getDistanceMatrix());
            distanceAverage+= currentDistance;
        }
        ancestralAverage = ancestralAverage/testData.size();
        distanceAverage = distanceAverage/testData.size();
        int n = mutationList.size();
        //upper bound of difference in distance is mutations all in one cluster compared to a chain (w/ one mutation per vertices)
        //this distance follows the patterns of tetrahedral numbers
        double tetrahedral = (n*(n+1)*(n+2))/6;
        double normalizedDistance = 1-(distanceAverage / tetrahedral);
        //creates probability between 0 and 1
        double probability = alpha*ancestralAverage + (1-alpha)*normalizedDistance;
        //System.out.println("probability: " + probability + " distance: " + normalizedDistance + " (" + distanceAverage + ") ancestral: " + ancestralAverage);
        return probability;
    }
    
    //calculates the beta probability for a tree based on input tree(s)
    //the probability is based on the probability that the population frequency of a cluster
    //is greater than the sum of the population frequencies of its children
    public static double betaProbability(Tree currentTree, ArrayList<String> mutationList, ArrayList<Tree> testTrees) {
        double probability = 1.0;
        double [] averageFrequencies = new double[mutationList.size()];
        for(Tree t:testTrees) {
            double[] curFrequencies = t.getFrequencies();
            for(int i=0; i<curFrequencies.length; i++) {
                averageFrequencies[i] = averageFrequencies[i] + curFrequencies[i];
            }
        }
        for(int i=0; i<averageFrequencies.length; i++) {
            averageFrequencies[i] = averageFrequencies[i]/testTrees.size();
        }
            int coverage = 100;
            for(int i=0; i<currentTree.getAdjacencyMatrix().length; i++) {
                double averageVariant = 0.0;
                double averageReference = 0.0;
                boolean hasChildren = false;
                for(String m:currentTree.getClusters().get(i)) {
                    int mutation1 = mutationList.indexOf(m);
                    double frequency = averageFrequencies[mutation1];
                    double variant = (frequency/100.0)*coverage;
                    double reference = (double) coverage - variant;
                    if(reference == 0) {
                        reference = .0001;
                    }
                    averageVariant += variant;
                    averageReference += reference;
                }
                averageVariant = averageVariant/currentTree.getClusters().get(i).size();
                averageReference = averageReference/currentTree.getClusters().get(i).size();
                double variantSum = 0.0;
                double referenceSum = 0.0;
                for(int j=0; j<currentTree.getAdjacencyMatrix().length; j++) {    
                    if(currentTree.getAdjacencyMatrix()[i][j]) {
                        hasChildren = true;
                        double averageVariant2 = 0.0;
                        double averageReference2 = 0.0;    
                        for(String m2:currentTree.getClusters().get(j)) {
                            int mutation2 = mutationList.indexOf(m2);
                            double frequency2 = averageFrequencies[mutation2];
                            double variant2 = (frequency2/100.0)*coverage;
                            double reference2 = (double) coverage - variant2;
                            if(reference2 == 0) {
                                reference2 = .0001;
                            }
                            averageVariant2 += variant2;
                            averageReference2 += reference2;
                        }
                        averageVariant2 = averageVariant2/currentTree.getClusters().get(j).size();
                        averageReference2 = averageReference2/currentTree.getClusters().get(j).size();
                        variantSum += averageVariant2;
                        referenceSum += averageReference2;
                    }
                }
                if(hasChildren) {
                    double comparison = betaComparison(averageVariant, averageReference, variantSum, referenceSum);
                    probability = probability * comparison;
                }

            }
            return probability;
    }

    /* implementation of Cook Beta Inequality algorithm for comparison
     of 2 beta distributed independent random variables
     based on code from AncesTree w/ help from Layla Oesper and Michael Hoffert
    */
    
    //h term for recurrence
    public static double h(double a, double b, double c, double d) {
        double hValue = 2;
        //Math.exp(org.apache.commons.math3.special.Beta.logBeta(a + c, b + d) - (org.apache.commons.math3.special.Beta.logBeta(a, b) + org.apache.commons.math3.special.Beta.logBeta(c, d)));
        return hValue;
        }

    //compare function: takes a, b, c, d values and outputs posterior probabilities
    public static double betaComparison(double a, double b, double c, double d) {
        //a = var mut1 
        //b = ref mut1
        //c = var mut2
        //d =  ref mut2
        double[] acList = {a, c};
        double[] bdList = {b, d};

        double aa = Math.min(acList[0], acList[1]);
        double bb = Math.min(bdList[0], bdList[1]);
        double cc = aa;
        double dd = bb;

        double result = 0.5;

        //iterating through values until g(a, b, c, d) is reached

        while (aa < acList[0]) {

            result += ( h(aa, bb, cc, dd) / aa);
            aa += 1;
        }

        while (bb < bdList[0]){

            result -= ( h(aa, bb, cc, dd) / bb);
            bb += 1;
            }
        while (cc < acList[1]){

            result -= ( h(aa, bb, cc, dd) / cc);
            cc += 1;
            } 
        while (dd < bdList[1]){

            result += ( h(aa, bb, cc, dd) / dd);
            dd += 1;
            }
        return result;

        }

    
    //find the mutations present in all of the data trees read in from a file
     public static ArrayList<String> commonMutations(String fileName) {
        ArrayList<String> mutationList = new ArrayList<String>();
        ArrayList<ArrayList<String>> allMutations = new ArrayList<ArrayList<String>>();
        try {
            File readFile = new File(fileName);
            Scanner sc = new Scanner(readFile);
            String curLine = sc.nextLine();
            while(!curLine.equals("#end")) {
                if(curLine.equals("#tree clusters")) {
                    ArrayList<String> currentMutations = new ArrayList<String>();
                    curLine = sc.nextLine();
                    while(curLine.length()>0 && curLine.charAt(0)!='#') {
                        String[] row = curLine.split("; ");
                        for(String m:row) {
                            currentMutations.add(m);
                        }
                        curLine = sc.nextLine();
                    }
                    allMutations.add(currentMutations);
                }
                curLine = sc.nextLine();
            }
            for(String m: allMutations.get(0)) {
                boolean inAll = true;
                for(int i=1; i<allMutations.size(); i++) {
                    if(!allMutations.get(i).contains(m)) {
                        inAll = false;
                    }
                }
                if(inAll) {
                    mutationList.add(m);
                }
            }
        }  catch(Exception ex) {
                    ex.printStackTrace();
        }
        return mutationList;
    }
    
    //creates trees for each of the sample trees read in from a file
    public static ArrayList<Tree> readAllDataFile(String fileName, ArrayList<String> mutationList) {
        ArrayList<Tree> testTrees = new ArrayList<Tree>();
        try {
            File readFile = new File(fileName);
            Scanner sc = new Scanner(readFile);
            String curLine = sc.nextLine();
            while(!curLine.equals("#end")){
                curLine = sc.nextLine();
                String[] number = curLine.split(" ");
                boolean[][] adjacencyMatrix = new boolean[number.length][number.length];
                ArrayList<ArrayList<String>> clusters = new ArrayList<ArrayList<String>>();
                int i =0;
                while(curLine.length()>0 && curLine.charAt(0)!='#') {
                        String[] row = curLine.split(" ");
                        for(int j=0; j<number.length; j++) {
                            if(Integer.parseInt(row[j])==1) {
                                adjacencyMatrix[i][j] = true;
                            }
                        }
                        curLine = sc.nextLine();
                        i++;
                    }
                while(curLine.length()==0) { 
                    curLine = sc.nextLine();
                }
                curLine = sc.nextLine();
                while(curLine.length()>0) {
                    ArrayList<String> group = new ArrayList<String>();
                    String[] parts = curLine.split("; ");
                    for(String c:parts) {
                        group.add(c);
                    }
                    clusters.add(group);
                    curLine=sc.nextLine();
                }
                curLine=sc.nextLine();
                curLine=sc.nextLine();
                double[] frequencies = new double[clusters.size()];
                String[] f = curLine.split(" ");
                for(int y=0; y<f.length; y++) {
                    frequencies[y] = Double.parseDouble(f[y])*100;
                }
                double[][] distanceMatrix = new double[1][1];
                Tree t = new Tree(clusters, adjacencyMatrix, clusters.size(), frequencies, distanceMatrix);
                testTrees.add(t);
                while(curLine.length()==0 || curLine.charAt(0)!='#') { 
                    curLine = sc.nextLine();
                }
            }
        }  catch(Exception ex) {
                    ex.printStackTrace();
        }
        return testTrees;
    }

    //creates trees containing only mutations found in all sample trees
    public static ArrayList<Tree> reduceTrees(ArrayList<Tree> testTrees, ArrayList<String> mutationList) {
        ArrayList<Tree> newTrees = new ArrayList<Tree>();
        for(Tree t: testTrees) {
            ArrayList<Integer> present = new ArrayList<Integer>();
            ArrayList<ArrayList<String>> newClusters = new ArrayList<ArrayList<String>>();
            double[] newFrequencies = new double[mutationList.size()];
            boolean[][] adjacencyMatrix = t.getAdjacencyMatrix();
            double[] frequencies = t.getFrequencies();
            ArrayList<ArrayList<String>> clusters = t.getClusters();
            for(int i =0; i<clusters.size(); i++) {
                ArrayList<String> group = new ArrayList<String>();
                for(String m:clusters.get(i)) {
                    if(mutationList.contains(m)) {
                        group.add(m);
                    }
                }
                if(group.size() !=0) {
                    newClusters.add(group);
                    present.add(i);
                }
            }
            boolean[][] newAdjacencyMatrix = new boolean[newClusters.size()][newClusters.size()];
           
            for(int j: present) {
                if(j!=0) {
                    for(int i=0; i<adjacencyMatrix.length; i++) {
                        if(adjacencyMatrix[i][j] && present.contains(i)) {
                            newAdjacencyMatrix[present.indexOf(i)][present.indexOf(j)] = true;
                        }
                        if(adjacencyMatrix[i][j] && !present.contains(i)) {
                            int parent = i;
                            while(!present.contains(parent)) {
                                for(int x=0; x<adjacencyMatrix.length; x++) {
                                    if(adjacencyMatrix[x][parent]) {
                                        parent = x;
                                    }
                                }
                            }
                            newAdjacencyMatrix[present.indexOf(parent)][present.indexOf(j)] = true;
                        }
                    }
                }
            }

            for(int i=0; i<mutationList.size(); i++) {
                int clusterNum = 0;
                for(int j=0; j<clusters.size(); j++) {
                    if(clusters.get(j).contains(mutationList.get(i))) {
                        clusterNum = j;
                    }
                }
                newFrequencies[i] = frequencies[clusterNum];
            }
            double[][] distanceMatrix = makeDistanceMatrix(mutationList, newClusters, newAdjacencyMatrix);
            Tree newTree = new Tree(newClusters, newAdjacencyMatrix, newClusters.size(), newFrequencies, distanceMatrix);
            newTrees.add(newTree);
        }
        return newTrees;
    }
    
    //DBMC can return a tree not rooted at the normal cell
    //change the direction of edges to make the tree rooted at the normal cell
    //if the normal cell is clustered with other mutations, make those mutations
    //the child of the normal cell and all other clusters descendants of that cluster
    public static Tree reroot(Tree t) {
        int root =0;
        int oldParent = 0;
        ArrayList<ArrayList<String>> cluster = new ArrayList<ArrayList<String>>();
        cluster = t.getClusters();
        boolean[][] aM = t.getAdjacencyMatrix();
        boolean[][] newAdjacencyMatrix = new boolean[aM.length][aM.length];
        for(int i=0; i<t.getClusters().size(); i++) {
            if(t.getClusters().get(i).contains("Normal Cell")) {
                root = i;
            }
        }
        for(int x=0; x<aM.length; x++) {
            for(int y=0; y<aM.length; y++) {
                if(y==root) {
                    if(aM[x][y]) {
                        newAdjacencyMatrix[y][x] = true;
                        oldParent = x;
                    }
                }
                else {
                    if(aM[x][y]) {
                        newAdjacencyMatrix[x][y] = true;
                    }
                }
            }
        }
        if(oldParent!=0) {
            boolean unchanged = false;
            while(!unchanged) {
                unchanged = true;
                for(int x=0; x<aM.length; x++) {
                    if(aM[x][oldParent]) {
                        newAdjacencyMatrix[oldParent][x] = true;
                        newAdjacencyMatrix[x][oldParent] = false;
                        oldParent=x;
                        unchanged = false;
                    }
                }
            }
        }
        //if the root contains mutations, make a new cluster containing those mutations 
        //and make that cluster the child of the root and all other cluster descendants
        //of the new cluster
        if(t.getClusters().get(root).size() > 1) {
            ArrayList<String> newCluster = new ArrayList<String>();
            for(int i=0; i<t.getClusters().get(root).size(); i++) {
                if(!t.getClusters().get(root).get(i).equals("Normal Cell")) {
                    newCluster.add(t.getClusters().get(root).get(i));
                }
            }
            cluster.add(newCluster);
            boolean [][] largerAdjacencyMatrix = new boolean[newAdjacencyMatrix.length + 1][newAdjacencyMatrix.length + 1];
        for(int i =0; i<newAdjacencyMatrix.length; i++) {
            for(int j=0; j<newAdjacencyMatrix.length; j++) {
                if(i==root) {
                    largerAdjacencyMatrix[newAdjacencyMatrix.length][j] = newAdjacencyMatrix[root][j] ;
                } else {
                largerAdjacencyMatrix[i][j] = newAdjacencyMatrix[i][j];
                }
            }
        }
            largerAdjacencyMatrix[root][newAdjacencyMatrix.length] = true;
            newAdjacencyMatrix = largerAdjacencyMatrix;
            
        }
        Tree newTree = new Tree(cluster, newAdjacencyMatrix, t.getClusters().size(), t.getFrequencies(), t.getDistanceMatrix());
        return newTree;
    }
    
    //calculates the rand index between two trees, which tells how similar the
    //clustering of the two trees is
    public static double randIndex(Tree tree1, ArrayList<Tree> treeList) {
        double rand = 1.0;
        for(Tree tree2:treeList) {
            double[][] distance1 = tree1.getDistanceMatrix();
            double[][] distance2 = tree2.getDistanceMatrix();
            double a = 0.0;
            double b = 0.0;
            double c = 0.0;
            double d = 0.0;
            for(int i=0; i<distance1.length-1; i++) {
                for(int j=i+1; j<distance1.length; j++) {
                    if(distance1[i][j] == 0 && distance2[i][j] == 0) {
                        a+=1.0;
                    }
                    else if(distance1[i][j] != 0 && distance2[i][j] != 0) {
                        b+=1.0;
                    }
                    else if(distance1[i][j] == 0 && distance2[i][j] != 0) {
                        c+=1.0;
                    } else {
                        d+=1.0;
                    }
                }
            }
            rand = rand * ((a+b)/(a+b+c+d));
        }
        return rand;
    }

    //creates a matrix a of all of the mutations where a[i][j] is true if 
    //mutation j is ancestral to mutation i
   public static boolean[][] ancestor(boolean [][] matrix, ArrayList<String> mutationList, ArrayList<ArrayList<String>> clusters) {
        ArrayList<Integer> parents = new ArrayList<Integer>();
        for(int i=0; i<matrix.length; i++) {
            boolean noParent = true;
            for(int j=0; j<matrix.length; j++) {
                if(matrix[j][i]){
                    noParent = false;
                }
            }
            if(noParent) {
                parents.add(i);
            }
        }
        int y =0;
        while(y < matrix.length && y<parents.size()) {
            for(int j =0; j<matrix.length; j++) {
                if(matrix[parents.get(y)][j]) {
                    parents.add(j);
                }
            }
            y++;
        }
       boolean [][] ancestors = new boolean[matrix.length][matrix.length];
       boolean[][] mutationAncestors = new boolean[mutationList.size()][mutationList.size()];
       for(int p:parents) {
           for(int j=0; j<matrix.length; j++) {
                if(matrix[j][p]) {
                    ancestors[p][j] = true;
                    for(int a=0; a<matrix.length; a++) {
                        if(ancestors[j][a]) {
                            ancestors[p][a] = true;
                        }
                    }
                }
            }
       }
      for(int c=0; c<clusters.size(); c++) {
        for(String m: clusters.get(c)) {
            int index = mutationList.indexOf(m);
            for(int z=0; z<ancestors.length; z++) {
                if(ancestors[c][z]) {
                    for(String m2: clusters.get(z)) {
                        int index2 = mutationList.indexOf(m2);
                        mutationAncestors[index][index2] = true;  
                    }
                }
               }
        }
      }
       
       return mutationAncestors;
    }
    
    //calculates the number of pairs of mutations m_1 and m_2 where m_2 is ancestral
    //to m_1 in tree 1 and in tree 2
    public static double ancestralDifferences(boolean[][] a1, boolean[][] a2) {
        double overlap = 0.0;
        double total = 0.0;
        int size = a1.length;
        if(a2.length<size){
            size = a2.length;
        }
        for(int i=0; i<size; i++) {
           for(int j=0; j<size; j++) {
                if(a1[i][j] && a2[i][j]) {
                    overlap +=1.0;
                } if(a1[i][j]) {
                    total +=1.0;
                }
            }
        }
        return overlap/total;
    }
    
    //creates a dotFile of a tree
    //outlines clusters in red if the sample frequencies indicate that it is less likely
    //that that cluster meets to the sum rule
    public static void writeDOT(boolean[][] writeAdjacencyMatrix, String fileName, ArrayList<ArrayList<String>> clusterNames, ArrayList<Integer> problems) {
        try {
			FileWriter f = new FileWriter(fileName);
			String newLine = System.getProperty("line.separator");
			f.write("digraph mytree {" + newLine);
			for (int i=0; i<writeAdjacencyMatrix.length; i++) {
                String parent = i + ": ";
                for(String name:clusterNames.get(i)) {
                    parent = parent + " " + name;
                }
                parent = "\"" + parent + "\"";
                boolean noChildren = true;
                for(int j=0; j<writeAdjacencyMatrix.length; j++) {
                    if(writeAdjacencyMatrix[i][j]) {
                        noChildren = false;
                        String child = j + ": ";
                        for(String childName:clusterNames.get(j)) {
                            child = child + " " + childName;
                        }
                        child = "\"" + child + "\"";
                        String dotLine = parent + " -> " + child + ";";
                        f.write(dotLine + newLine);
                        if(problems.contains(i)) {
                            String conflict = parent + " [color=red];";
                            f.write(conflict + newLine);
                        }
                    }
                }
                if(noChildren) {
                    parent = parent + ";";
                    f.write(parent + newLine);
                }

            }			
			f.write("}" + newLine);
			f.close();
		}
		catch(Exception ex) {
				ex.printStackTrace();
		}
	}
    
    //outputs the complete-data distance, distance to true tree, ancestral similarity with true tree,
    //clustering similarity with true tree, and beta probabilities for all trials with simulated data
    //writes it to a file called outputFile.txt
    public static void printEvaluations(double [][] matrix) {
         try {
			FileWriter f = new FileWriter("dotTest.csv");
			String newLine = System.getProperty("line.separator");
            f.write("Output Data" + newLine);
            String explanation = "True Tree Distance Sum, MV Distance Sum, DBMC Distance Sum, MCMC Distance Sum, MV Distance to, , DBMC Distance to, , MCMC distance To, , ";
             explanation = explanation + "sample 1 distance, sample 2 distance, sample 3 distance, sample 4 distance, sample 5 distance, ";
             explanation = explanation + "MV ancestral similarity, , DBMC ancestral similarity, ,  MCMC ancestral similarity, , sample 1 ancestral, sample 2 ancestral, sample 3 ancestral, sample 4 ancestral, sample 5 ancestral, ";
             explanation = explanation + "MV cluster similarity, , DBMC cluster similarity, , MCMC cluster similarity, , sample 1 clustering, sample 2 clustering, sample 3 clustering, sample 4 clustering, sample 5 clustering, ";
             explanation = explanation + "MV Beta probability, DBMC Beta probability, MCMC Beta Probability, ";
             explanation = explanation + "MV New Metric (True), , DBMC New Metric (True), , MCMC New Metric (true), , sample 1 new, sample 2 new, sample 3 new, sample 4 new, sample 5 new, ";
             explanation = explanation + "MV New Metric, DBMC New Metric, MCMC New Metric";
             
             f.write(newLine);
             for(int i=0; i<matrix.length; i++) {
                String line = "";
                 for(int j=0; j<matrix[0].length; j++) {
                     line = line + matrix[i][j] + ", ";
                    }
                 f.write(line + newLine);
             }
            f.close();
            }
            catch(Exception ex) {
				ex.printStackTrace();
		}
    }
    //creates matrix containing all of the output information for each trial
    public static double[][] outputMatrix(Tree trueTree, Tree majorityVoteTree, Tree mcDistanceTree, ArrayList<Tree> testData, ArrayList<String> mutationList, int trial, double [][] matrix, Tree mcTree1, Tree mcTree2, Tree mcTree3, Tree mcTree4, Tree mcTree5) {
            ArrayList<Tree> trueList = new ArrayList<Tree>();
            trueList.add(trueTree);
             boolean[][] a1 = ancestor(trueTree.getAdjacencyMatrix(), mutationList, trueTree.getClusters());
            double ancestralAverage = 0.0;
            double randAverage = 0.0;
            double distanceAverage = 0.0;
            int counter = 0;
            double al =0;
            for(Tree t:testData) {
                boolean [][] an = ancestor(t.getAdjacencyMatrix(), mutationList, t.getClusters());
                double currentAncestral = ancestralDifferences(a1, an);
                ancestralAverage += currentAncestral;
                double currentRand = randIndex(t, trueList);
                randAverage += currentRand;
                int currentDistance = calculateDistance(trueTree.getDistanceMatrix(), t.getDistanceMatrix());
                double currentCombined1 = combinedComparison(trueList, t, mutationList, 0.0);
                double currentCombined2 = combinedComparison(trueList, t, mutationList, 0.3);
                double currentCombined3 = combinedComparison(trueList, t, mutationList, 0.5);
                double currentCombined4 = combinedComparison(trueList, t, mutationList, 0.7);
                double currentCombined5 = combinedComparison(trueList, t, mutationList, 1.0);
                if(counter <=4) {
                    matrix[trial][15+counter] = currentDistance; 
                    matrix[trial][27+counter] = currentAncestral; 
                    matrix[trial][39+counter] = currentRand;
                    matrix[trial][51+counter] = currentCombined1;
                    matrix[trial][63+counter] = currentCombined2;
                    matrix[trial][75+counter] = currentCombined3;
                    matrix[trial][87+counter] = currentCombined4;
                    matrix[trial][99+counter] = currentCombined5;

                }
                distanceAverage += currentDistance;
                counter++;
            }
            boolean [][] a2 = ancestor(majorityVoteTree.getAdjacencyMatrix(), mutationList, majorityVoteTree.getClusters());
            boolean [][] a3 = ancestor(mcDistanceTree.getAdjacencyMatrix(), mutationList, mcDistanceTree.getClusters());
            boolean [][] a4 = ancestor(mcTree1.getAdjacencyMatrix(), mutationList, mcTree1.getClusters());
            boolean [][] a5 = ancestor(mcTree2.getAdjacencyMatrix(), mutationList, mcTree2.getClusters());
            boolean [][] a6 = ancestor(mcTree3.getAdjacencyMatrix(), mutationList, mcTree3.getClusters());
            boolean [][] a7 = ancestor(mcTree4.getAdjacencyMatrix(), mutationList, mcTree4.getClusters());
            boolean [][] a8 = ancestor(mcTree5.getAdjacencyMatrix(), mutationList, mcTree5.getClusters());

           
			matrix[trial][0] = distanceSum(testData, trueTree);
            matrix[trial][1] = distanceSum(testData, majorityVoteTree);
            matrix[trial][2] = distanceSum(testData, mcDistanceTree);
            matrix[trial][3] = distanceSum(testData, mcTree1);
            matrix[trial][4] = distanceSum(testData, mcTree2);
            matrix[trial][5] = distanceSum(testData, mcTree3);
            matrix[trial][6] = distanceSum(testData, mcTree4);
            matrix[trial][7] = distanceSum(testData, mcTree5);
            matrix[trial][8] = calculateDistance(trueTree.getDistanceMatrix(), majorityVoteTree.getDistanceMatrix());
            matrix[trial][9] = calculateDistance(trueTree.getDistanceMatrix(), mcDistanceTree.getDistanceMatrix());
            matrix[trial][10] = calculateDistance(trueTree.getDistanceMatrix(), mcTree1.getDistanceMatrix());
            matrix[trial][11] = calculateDistance(trueTree.getDistanceMatrix(), mcTree2.getDistanceMatrix());
            matrix[trial][12] = calculateDistance(trueTree.getDistanceMatrix(), mcTree3.getDistanceMatrix());
            matrix[trial][13] = calculateDistance(trueTree.getDistanceMatrix(), mcTree4.getDistanceMatrix());
            matrix[trial][14] = calculateDistance(trueTree.getDistanceMatrix(), mcTree5.getDistanceMatrix());
            matrix[trial][20] = ancestralDifferences(a1, a2);
            matrix[trial][21] = ancestralDifferences(a1, a3);
            matrix[trial][22] = ancestralDifferences(a1, a4);
            matrix[trial][23] = ancestralDifferences(a1, a5);
            matrix[trial][24] = ancestralDifferences(a1, a6);
            matrix[trial][25] = ancestralDifferences(a1, a7);
            matrix[trial][26] = ancestralDifferences(a1, a8);
            matrix[trial][32] = randIndex(majorityVoteTree, trueList);
            matrix[trial][33] = randIndex(mcDistanceTree, trueList);
            matrix[trial][34] = randIndex(mcTree1, trueList);
            matrix[trial][35] = randIndex(mcTree2, trueList);
            matrix[trial][36] = randIndex(mcTree3, trueList);
            matrix[trial][37] = randIndex(mcTree4, trueList);
            matrix[trial][38] = randIndex(mcTree5, trueList);

            matrix[trial][44] = combinedComparison(trueList, majorityVoteTree, mutationList, 0);
            matrix[trial][45] = combinedComparison(trueList, mcDistanceTree, mutationList, 0);
            matrix[trial][46] = combinedComparison(trueList, mcTree1, mutationList, 0);
            matrix[trial][47] = combinedComparison(trueList, mcTree2, mutationList, 0);
            matrix[trial][48] = combinedComparison(trueList, mcTree3, mutationList, 0);
            matrix[trial][49] = combinedComparison(trueList, mcTree4, mutationList, 0);
            matrix[trial][50] = combinedComparison(trueList, mcTree5, mutationList, 0);

            matrix[trial][56] = combinedComparison(trueList, majorityVoteTree, mutationList, .3);
            matrix[trial][57] = combinedComparison(trueList, mcDistanceTree, mutationList, .3);
            matrix[trial][58] = combinedComparison(trueList, mcTree1, mutationList, .3);
            matrix[trial][59] = combinedComparison(trueList, mcTree2, mutationList, .3);
            matrix[trial][60] = combinedComparison(trueList, mcTree3, mutationList, .3);
            matrix[trial][61] = combinedComparison(trueList, mcTree4, mutationList, .3);
            matrix[trial][62] = combinedComparison(trueList, mcTree5, mutationList, .3);

            matrix[trial][68] = combinedComparison(trueList, majorityVoteTree, mutationList, .5);
            matrix[trial][69] = combinedComparison(trueList, mcDistanceTree, mutationList, .5);
            matrix[trial][70] = combinedComparison(trueList, mcTree1, mutationList, .5);
            matrix[trial][71] = combinedComparison(trueList, mcTree2, mutationList, .5);
            matrix[trial][72] = combinedComparison(trueList, mcTree3, mutationList, .5);
            matrix[trial][73] = combinedComparison(trueList, mcTree4, mutationList, .5);
            matrix[trial][74] = combinedComparison(trueList, mcTree5, mutationList, .5);

            matrix[trial][80] = combinedComparison(trueList, majorityVoteTree, mutationList, .7);
            matrix[trial][81] = combinedComparison(trueList, mcDistanceTree, mutationList, .7);
            matrix[trial][82] = combinedComparison(trueList, mcTree1, mutationList, .7);
            matrix[trial][83] = combinedComparison(trueList, mcTree2, mutationList, .7);
            matrix[trial][84] = combinedComparison(trueList, mcTree3, mutationList, .7);
            matrix[trial][85] = combinedComparison(trueList, mcTree4, mutationList, .7);
            matrix[trial][86] = combinedComparison(trueList, mcTree5, mutationList, .7);

            matrix[trial][92] = combinedComparison(trueList, majorityVoteTree, mutationList, 1);
            matrix[trial][93] = combinedComparison(trueList, mcDistanceTree, mutationList, 1);
            matrix[trial][94] = combinedComparison(trueList, mcTree1, mutationList, 1);
            matrix[trial][95] = combinedComparison(trueList, mcTree2, mutationList, 1);
            matrix[trial][96] = combinedComparison(trueList, mcTree3, mutationList, 1);
            matrix[trial][97] = combinedComparison(trueList, mcTree4, mutationList, 1);
            matrix[trial][98] = combinedComparison(trueList, mcTree5, mutationList, 1);
        
        return matrix;
    }
}