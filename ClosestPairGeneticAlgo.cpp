/*
 * Luke Duncan
 *
 * CIS 405 - Term Project
 *
 * Description:
 * Implementation and comparision of a genetic closest pair of points algorithm vs. a divide and conquer
 * algorithm
 *
 * Important Note:
 * While 100% of this code was written by myself, I have submitted the divide-and-conquer closest pair
 * algorithm as an assignment from another class this semester (CIS 306 - Discrete II).  This project
 * is an extension of that assignment.  I took my interest in AI algorithms and decided to implement
 * a genetic algo solution for a problem discussed in class.
 */

// HEADERS /////////////////////////////////////////////////////////
#include<iostream>
#include<vector>
#include<algorithm>
#include<cmath>
#include<cstdio>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<limits>
#include<set>
using namespace std;

#define POPULATIONSIZE 10
#define MAXLOOP 1000
#define CROSSOVERMAX 6
#define NUMOFCROSSOVERS 2
#define POTENTIALMUTATIONS 2
#define NUMOFRUNS 40

// DATA STRUCTURES /////////////////////////////////////////////////

struct point
{
    double x, y;
    int index;
};

struct chromosome
{
    point points[2];
    double fitness;
};

double calcDistance(point p1, point p2);
void splitVector(vector<point> allPoints, vector<point> &left, vector<point> &right, int line);
double closestPairPlaneSweep(vector<point> allPoints, double delta);
double closestPair(vector<point> allPoints);

bool sortByFitnessDesc (chromosome chrome1, chromosome chrome2) { return (chrome1.fitness > chrome2.fitness); }
bool sortByFitnessAsc (chromosome chrome1, chromosome chrome2) { return (chrome1.fitness < chrome2.fitness); }

// CUSTOM SORTS ////////////////////////////////////////////////////
bool sortByX (point p1, point p2) { return (p1.x < p2.x); }
bool sortByY (point p1, point p2) { return (p1.y < p2.y); }
bool sortByXY (point p1, point p2)
{
    if(p1.x == p2.x) return (p1.y < p2.y);
    else return (p1.x < p2.x);
}

// CLOSEST PAIR ALGO'S /////////////////////////////////////////////
double calcDistance(point p1, point p2)
/*
 * Pre:        Two points are passed
 * Post:    The distance between them is calculated using Pathagorean Theorem
 */
{
    double a;
    double b;

    a = p1.x - p2.x;
    b = p1.y - p2.y;

    return sqrt(a*a + b*b);
}

void splitVector(vector<point> allPoints, vector<point> &left, vector<point> &right, int line)
/*
 * Pre:        A vector X is passed, with two empy vectors L and R
 * Post:    The left half of X is stored in L and the right in R.
 */
{
    for(int i=0; i<allPoints.size(); i++)
    {
        if(i<line) left.push_back(allPoints[i]);
        else right.push_back(allPoints[i]);
    }

}

double closestPairPlaneSweep(vector<point> allPoints, double delta)
/*
 * Pre:        A vector of points, sorted by Y-coordinates, is passed by Divide-and-Conquer
 *             algo
 * Post:    The shortest distance between any two points is returned as a double
 */
{
    for(int i=0; i<allPoints.size(); i++)
    {
        for(int j=i+1; j<allPoints.size() && (allPoints[j].y > (allPoints[i].y - delta)); j++)
        // Check everything in front of the current point
        // stoping when the Sweeping line is greater than delta away from the point
        {
            double curDelta = calcDistance(allPoints[i], allPoints[j]);
            delta = min(delta, curDelta);
        }
    }
    return delta;
}

double closestPair(vector<point> allPoints)
/*
 * Pre:        A vector of points, sorted by X-Coordinates, is passed
 * Post:    A Divide-and-Conquer closest pair algorithm calculates the
 *             the smallest distance between any two points.  A double is returned
 */
{
    if(allPoints.size() == 1) return 10000;

    int seperationLine = (allPoints.size())/2;
    vector<point> left;
    vector<point> right;
    double deltaLeft;
    double deltaRight;
    double delta;
    double deltaMid;
    double deltaMidLeft;
    double deltaMidRight;

    // Split the vector in half, recursively find
    // closests pair
    splitVector(allPoints, left, right, seperationLine);
    deltaLeft = closestPair(left);
    deltaRight = closestPair(right);
    delta = min(deltaLeft, deltaRight);    // Store min for use in mid section

    // Calculate Mid Section
    // for plane-sweep algo
    deltaMid = double(allPoints[seperationLine].x);
    deltaMidLeft = deltaMid - delta;
    deltaMidRight = deltaMid + delta;

    while(!(allPoints[0].x >= deltaMidLeft && allPoints[0].x <= deltaMidRight))
    // Delete everything to the left and not in delta
    {
        allPoints.erase(allPoints.begin());
    }

    while(!(allPoints[allPoints.size()-1].x >= deltaMidLeft && allPoints[allPoints.size()-1].x <= deltaMidRight))
    // Delete everything to the right and not in delta
    {
        allPoints.erase(allPoints.end()-1);
    }

    // Sort by y for a y-directional plane sweep
    sort(allPoints.begin(), allPoints.end(), sortByY);

    // return result of plane-sweep
    return closestPairPlaneSweep(allPoints, delta);
}

// CLOSEST PAIR GENETIC ALGO ///////////////////////////////////////

void printPoint(point passed)
/*
* Pre:    a point is passed
* Post: it's co-ordinates are printed.  Used for debugging
*/
{
    cout << "("<< passed.x << ", " << passed.y<<")\n";
}

void printPopulation(vector<chromosome> &population)
/*
* Pre:    a population is passed
* Post: The chromosomes are printed.  Used for debugging
*/
{
    cout << "Population Status:\n";
    for(int i=0; i<population.size(); i++)
    {
        point p1, p2;
        double p1x, p1y, p2x, p2y;
        p1 = population[i].points[0];
        p2 = population[i].points[1];
        cout << "Chromosome " << i << ":\n";
        cout << "{(" << p1.x << "," << p1.y << "),(" << p2.x<<", " << p2.y<<")} Fitness: " << population[i].fitness << endl;
        cout << "Indices: " << p1.index << ", " << p2.index << endl;
    }

    cout << endl <<  "---------------------------------------------------\n\n";
}

double fitnessPopulation(vector<chromosome> &population)
/*
* Pre:    a population is passed
* Post: The fitness of each chromosome is determined
*/
{
    double minDistance = numeric_limits<double>::infinity();
    vector<double> popFitness;
    for(int i=0; i< population.size(); i++)
    {
        double cur = calcDistance(population[i].points[0], population[i].points[1]);
        population[i].fitness = cur;
        if(cur < minDistance) minDistance = cur;
    }

    return minDistance;
}

double fitnessChromosome(chromosome &chrome1)
/*
* Pre:    a chromosome is passed
* Post: it's fitness is determined
*/
{
    chrome1.fitness = calcDistance(chrome1.points[0],chrome1.points[0]);
}

void mutation(chromosome &chrome1, vector<point> &allPoints)
/*
* Pre:    a chromosome, and all possible points are passed
* Post: a random mutation is executed
*/
{
    int toReplace, toGet;
    toReplace = rand() % 2; // 0 or 1, x or y
    toGet = rand() % allPoints.size();

    if(toReplace  && chrome1.points[1].index != allPoints[toGet].index) chrome1.points[0] = allPoints[toGet];
    else if (!toReplace  && chrome1.points[0].index != allPoints[toGet].index) chrome1.points[1] = allPoints[toGet];
}

void calcMutation(vector<chromosome> &population, vector<chromosome> &nextGen, vector<point> &allPoints)
/*
* Pre:    a population is passed
* Post: a chromosome is selected for mutation and mutated
*/
{
    int toReplace = rand() % population.size();
    chromosome chrome1 = population[toReplace];
    mutation(chrome1, allPoints);
    nextGen.push_back(chrome1);
}

void crossOver(chromosome &chrome1, chromosome &chrome2)
/*
* Pre:    two chromosomes are passed
* Post: a crossover is produced
*/
{
    int toReplace, toGet;
    int toMutate = rand() % 2;
    toReplace = rand() % 2; // 0 or 1, x or y
    toGet = rand() % 2; // 0 or 1, x or y

    while(chrome1.points[toReplace].index == chrome1.points[toGet].index)
    {
        toReplace = rand() % 2; // 0 or 1, x or y
        toGet = rand() % 2; // 0 or 1, x or y
    }

    swap(chrome1.points[toReplace], chrome2.points[toGet]);
}

void calcCrossOver(vector<chromosome> &population, vector<chromosome> nextGen, vector<point> allPoints)
/*
* Pre:    a population is passed
* Post: chromosomes are selected for crossover and crossed.  Mutations can occur
*/
{
    // May want to allow it to pick two of the top four?
    int toMutate = rand() % 100;
    int toGet = rand() % CROSSOVERMAX;
    int toReplace = rand() % CROSSOVERMAX;
    while(toReplace == toGet) toReplace = rand() % CROSSOVERMAX;


    chromosome chrome1;
    chromosome chrome2;

    sort(population.begin(), population.end(), sortByFitnessDesc);
    chrome1 = population[toGet];
    chrome2 = population[toReplace];

    if(toMutate < 91) calcMutation(population, nextGen, allPoints);
    crossOver(chrome1, chrome2);
    nextGen.push_back(chrome1);
    nextGen.push_back(chrome2);
}

void generatePopulation(vector<chromosome> &population, vector<point> &allPoints)
/*
* Pre:    a population is passed
* Post: cross-overs and mutations occur to create the next generation
*/
{
    vector<chromosome> nextGen;
    // Do cross-overs
    for(int i=0; i < NUMOFCROSSOVERS; i++) calcCrossOver(population, nextGen, allPoints);

    // introduce mutations
//    for(int i=0; i < POTENTIALMUTATIONS; i++)
    calcMutation(population, nextGen, allPoints);

    // add nextGen to population
    for(int i=0; i<nextGen.size(); i++)
    {
        population.push_back(nextGen[i]);
    }
    fitnessPopulation(population);
}

void initializePopulation(vector<chromosome> &population, vector<point> &allPoints)
/*
* Pre:    an empty population is passed
* Post: a random population is generated
*/
{
    while(population.size() != POPULATIONSIZE)
    {
        chromosome chromeCur;
        vector<point> curPoints;
        int previous = -1;

        int spot = 0;
        while(curPoints.size() != 2)
        {
            point pointCur;

            int toGet = rand() % allPoints.size();
            while(toGet == previous) toGet = rand() % allPoints.size(); // ensure uniqueness of point

            previous = toGet;
            pointCur = allPoints[toGet];
            //pointCur.index = toGet;
            curPoints.push_back(pointCur);
        }


        chromeCur.points[0] = curPoints[0];
        chromeCur.points[1] = curPoints[1];

        population.push_back(chromeCur);
    }
    fitnessPopulation(population);
}

void survivalOfTheFittest(vector<chromosome> &population)
/*
* Pre:    a population is passed
* Post: the fittest members of the population survive. i.e. Darwinian computer science!
*/
{
    int toKill = population.size() - POPULATIONSIZE;

    sort(population.begin(), population.end(), sortByFitnessAsc);

    for(int i=0; i< toKill; i++)
    {
        population.erase(population.begin());
    }
}

double closestPairGenetic(vector<point> allPoints)
/*
* Pre:    a set of points is passed
* Post: a probabalistically correct shortest distance between points is returned.  Could
    also be considered optimized.
*/
{
    srand ( time(NULL) );
    vector<chromosome> population;
    double minDistance = numeric_limits<double>::infinity();
    int countSinceChange = -1;

    // initialize
    initializePopulation(population, allPoints);

//    cout << "Initial Population:\n";
//    printPopulation(population);

    while(countSinceChange != MAXLOOP)
    // loop until fitness doesn't change for X loops
    {
        countSinceChange++;

        // get fitness of population
        double populationMin = fitnessPopulation(population);
        if(populationMin < minDistance) {
            minDistance = populationMin;
            countSinceChange = 0;
        }

//        cout << "\n\nPop Vs. Overall: " << populationMin << "," << minDistance << "\n\n";

        // if population size is too big, kill the weakest chromosomes
        if(population.size() > POPULATIONSIZE)
        {
            survivalOfTheFittest(population);
        }

        // generation next population
        generatePopulation(population, allPoints);
//        printPopulation(population);
    }


    return minDistance;
}

// DRIVER AUX FUNCTIONS /////////////////////////////////////////

void printSet(vector<point> &allPoints)
/*
* Pre:    a set of points is passed
* Post: all of the points are printed
*/
{
    cout << "{";
    for(int i=0; i<allPoints.size(); i++)
    {
        cout << "("<< allPoints[i].x<< ", "<< allPoints[i].y<<"),";
    }
    cout << "}\n";
}

int getPercentage(int top, int bottom)
/*
* Pre:    a fraction is passed
* Post: an integer value is returned of it's percentage
*/
{
    double cur = double(top)/double(bottom);

    return int(cur*100);
}

// DRIVER ///////////////////////////////////////////////////////

int main()
{
    srand ( time(NULL) );

    for(int k=0; k<25; k++)
    {

        int total = 0;
        vector<point> allPoints;
        set<double> allReturned;
        double correct;
        int correctCount =0;
        double average =0;

        while(total < POPULATIONSIZE) total = rand() % 100;

        for(int i=0; i<total; i++)
        // Generate Random Points on a cartesian graph
        {
            point curPoint;
            curPoint.index = i; // necessary for genetic algo
            curPoint.x = rand() % 1000;
            curPoint.y = rand() % 1000;
            allPoints.push_back(curPoint);
        }

        cout << total << " random points were generated and are displayed below.  None of these points are garanteed to be unique.\n";
        printSet(allPoints);
        cout << endl << endl;


        sort(allPoints.begin(), allPoints.end(), sortByX);
        cout << "Correct shortest distance between points:\n";
        correct = closestPair(allPoints);
        cout << correct << endl << endl;

        for(int i=0; i<NUMOFRUNS; i++)
        {
            double cur = closestPairGenetic(allPoints);
            allReturned.insert(cur);
            if(cur == correct) correctCount++;
    //        cout << cur << endl;
            average += cur;
        }

        cout << "The genetic algorithm produced " << allReturned.size() << " distinct answer and was correct " << getPercentage(correctCount, NUMOFRUNS) <<"% of the time\n";
        cout << "The genetic algorithms average was off by " << (average/NUMOFRUNS) - correct << endl;

    }
    return 0;
}
