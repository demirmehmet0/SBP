#include <math.h>
#include <ctype.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include "random"
#include "array"
#include <sstream>

#define SIZE 16
#define TIME_LIMIT 90000
using namespace std;

const int MaxBaseSize = 1000;
const bool PRINTROWS = true;
struct Element
{
    int parent_i;
    int parent_j;
    int newDist[SIZE];
};
clock_t startProc;
int iter = 0;
int NumInputs;
int DepthLimit;
int NumMatrices;
int MaxDist;
int St;
long long int DepthNewBase;
int NumTargets;
int ProgramSize;
long long int Target[MaxBaseSize];
int Dist[MaxBaseSize];
int NDist[MaxBaseSize];
long long int Base[MaxBaseSize];
int BaseSize;
int TargetsFound;
char Result[MaxBaseSize][50];
int Res;
char *flag;
int Depth[MaxBaseSize];
int MaxDepth;
int minXOR = 999;
int chosenParam = 7;
long long int NewBase;

void InitBase();

void ReadTargetMatrix();

bool is_target(long long int x);

bool is_base(long long int x);

int NewDistance(int u);

int TotalDistance();

bool reachablE(long long int T, int K, int S, long long int L);

bool EasyMove();

void PickNewBaseElement();

void binprint(long long int x);

void PrintExpression(int No);

void PrintBase(void);

void PrintResult();

int Max(int a, int b);

string name = "Joltik";
ifstream TheMatrix;
mt19937 rand_generator;
default_random_engine generator;

int main(int argc, char *argv[])
{

    clock_t start = clock();
    while ((clock() - start) / CLOCKS_PER_SEC < TIME_LIMIT)
    {
        int Threshold;
        MaxDist = 0;
        DepthNewBase = 0;
        BaseSize = 0;
        TargetsFound = 0;
        MaxDepth = 0;
        int i = 0;

        clock_t t1 = clock();
        TheMatrix.open("" + name + ".txt");

        TheMatrix >> NumMatrices;
        TheMatrix >> Threshold;
        TheMatrix >> DepthLimit;

        for (i = 0; i < NumMatrices; i++)
        {
            startProc = clock();
            ReadTargetMatrix();
            if (MaxDist + 1 > pow(2, DepthLimit))
                continue;
            InitBase();
            ProgramSize = 0;
            int counter = 0;
            St = i + 1;

            while (TargetsFound < NumTargets)
            {

                counter++;
                if (ProgramSize > Threshold)
                {
                    break;
                }
                if (!EasyMove())
                    PickNewBaseElement();
            }

            if (ProgramSize <= Threshold)
                PrintResult();
            else
                printf("");
        }
        TheMatrix.close();
    }
}

void InitBase()
{
    TargetsFound = 0;
    Res = 0;
    Base[0] = 1;
    Depth[0] = 0;
    MaxDepth = 0;
    for (int i = 1; i < NumInputs; i++)
    {
        Base[i] = 2 * Base[i - 1];
        Depth[i] = 0;
    }
    BaseSize = NumInputs;
    for (int i = 0; i < NumTargets; i++)
        if (Dist[i] == 0)
        {
            TargetsFound++;

            for (int j = 0; j < NumInputs; ++j)
                if (Base[j] == Target[i])
                {
                    sprintf(Result[Res], "y%d = x%d  *  (0)\n", i, j);
                    ++Res;
                    break;
                }
        }
}

int TotalDistance()
{
    int D = 0;
    int t;
    for (int i = 0; i < NumTargets; i++)
    {
        t = NewDistance(i);
        NDist[i] = t;
        D = D + t;
    }
    return D;
}

bool EasyMove()
{
    int t;
    bool foundone = false;

    for (int i = 0; i < NumTargets; i++)
        if (Dist[i] == 1)
        {
            foundone = true;
            t = i;
            break;
        }
    if (!foundone)
        return false;

    NewBase = Target[t];
    Base[BaseSize] = NewBase;
    BaseSize++;
    DepthNewBase = pow(2, DepthLimit);
    for (int i = 0; i < BaseSize; ++i)
        for (int j = i + 1; j < BaseSize; ++j)
            if (((Base[i] ^ Base[j]) == Base[BaseSize - 1]) && (DepthNewBase > pow(2, Max(i, j) + 1)))
            {
                DepthNewBase = pow(2, Max(i, j) + 1);
            }

    for (int u = 0; u < NumTargets; u++)
        Dist[u] = NewDistance(u);
    ProgramSize++;

    PrintExpression(t);
    TargetsFound++;
    return true;
}

void PrintResult()
{
    std::fstream f("" + name + ".txt",
                   std::fstream::out | std::fstream::app);

    f << St << endl
      << endl;
    f << ProgramSize << endl
      << endl;
    if (ProgramSize <= minXOR)
    {
        minXOR = ProgramSize;
    }
    iter++;
    for (int i = 0; i < Res; ++i)
        f << Result[i];

    cout << "XOR Count= " << ProgramSize << "   Depth= " << MaxDepth << endl;
}

int Max(int a, int b)
{
    if (Depth[a] > Depth[b])
        return Depth[a];
    else
        return Depth[b];
}

void PrintExpression(int No)
{
    int i, j;

    for (i = 0; i < BaseSize; ++i){
        for (j = i + 1; j < BaseSize; ++j){
            if (((Base[i] ^ Base[j]) == Base[BaseSize - 1]) && (pow(2, Max(i, j) + 1) == DepthNewBase))
            {
                Depth[BaseSize - 1] = Max(i, j) + 1;
                if (Depth[BaseSize - 1] > MaxDepth)
                    MaxDepth = Depth[BaseSize - 1];
                flag = Result[Res];
                flag += sprintf(flag, "t%d = ", ProgramSize);
                if (i < NumInputs)
                    flag += sprintf(flag, "x%d + ", i);
                else
                    flag += sprintf(flag, "t%d + ", i - NumInputs + 1);
                if (j < NumInputs)
                    flag += sprintf(flag, "x%d *  y%d  (%d)\n", j, No, Depth[BaseSize - 1]);
                else
                    flag += sprintf(flag, "t%d *  y%d  (%d)\n", j - NumInputs + 1, No, Depth[BaseSize - 1]);
                ++Res;
                return;
            }
        }
    }
}

void PickNewBaseElement()
{
    int MinDistance;
    long long int TheBest;
    long long int TheBestArray[MaxBaseSize * NumTargets];
    int ThisDist;
    int ThisNorm, OldNorm;
    int besti, bestj, d;
    bool easytarget;
    int BestDist[MaxBaseSize];
    Element *AllElements = new Element[BaseSize * (BaseSize - 1)];
    int counter = 0;
    int count2 = 0;
    MinDistance = BaseSize * NumTargets;
    OldNorm = 0;

    for (int i = 0; i < BaseSize - 1; i++)
    {
        if (Depth[i] + 1 >= DepthLimit)
            continue;
        for (int j = i + 1; j < BaseSize; j++)
        {
            if (Depth[j] + 1 >= DepthLimit)
                continue;
            NewBase = Base[i] ^ Base[j];
            if (NewBase == 0)
            {
                cout << "a base is 0, should't happen " << endl;
                exit(0);
            }
            if (is_base(NewBase))
                continue;
            easytarget = false;
            if (is_target(NewBase))
            {
                cout << "shouldn't find an easy target here " << endl;
                exit(0);
            }
            DepthNewBase = pow(2, Max(j, i) + 1);
            ThisDist = TotalDistance();
            if (ThisDist <= MinDistance)
            {
                ThisNorm = 0;
                for (int k = 0; k < NumTargets; k++)
                {
                    d = NDist[k];
                    ThisNorm = ThisNorm + d * d;
                }
                if ((ThisDist < MinDistance) || (ThisNorm > OldNorm))
                {
                    if (counter > chosenParam)
                        counter = 0;
                    AllElements[counter].parent_i = i;
                    AllElements[counter].parent_j = j;
                    for (int uu = 0; uu < NumTargets; uu++)
                        AllElements[counter].newDist[uu] = NDist[uu];
                    TheBestArray[counter] = NewBase;
                    MinDistance = ThisDist;
                    OldNorm = ThisNorm;
                    counter++;
                    count2++;
                }
            }
        }
        if (easytarget)
            break;
    }

    rand_generator.seed(time(0));
    uniform_int_distribution<int> rand_distribution(0, counter - 1);
    int number = rand_distribution(rand_generator);
    if (counter == 0 || counter == 1)
        number = 0;
    if (counter == 2 && count2 == 2)
        number = 1;
    if (counter == 3 && count2 == 3)
        number = 2;
    besti = AllElements[number].parent_i;
    bestj = AllElements[number].parent_j;
    TheBest = TheBestArray[number];
    for (int uu = 0; uu < NumTargets; uu++)
        BestDist[uu] = AllElements[number].newDist[uu];

    NewBase = TheBest;
    for (int i = 0; i < NumTargets; i++)
        Dist[i] = BestDist[i];

    Base[BaseSize] = TheBest;
    Depth[BaseSize] = Max(besti, bestj) + 1;
    if (Depth[BaseSize] > MaxDepth)
        MaxDepth = Depth[BaseSize];
    BaseSize++;

    ProgramSize++;

    flag = Result[Res];
    flag += sprintf(flag, "t%d = ", ProgramSize);
    if (besti < NumInputs)
        flag += sprintf(flag, "x%d + ", besti);
    else
        flag += sprintf(flag, "t%d + ", besti - NumInputs + 1);
    if (bestj < NumInputs)
        flag += sprintf(flag, "x%d  (%d)\n", bestj, Depth[BaseSize - 1]);
    else
        flag += sprintf(flag, "t%d  (%d)\n", bestj - NumInputs + 1, Depth[BaseSize - 1]);
    ++Res;
    if (is_target(TheBest))
        TargetsFound++;
}

void binprint(long long int x)
{
    long long int t = x;
    for (int i = 0; i < NumInputs; i++)
    {
        if (t % 2)
            cout << "1 ";
        else
            cout << "0 ";
        t = t / 2;
    }
}

void PrintBase()
{
    int i;

    for (i = 0; i < BaseSize; ++i)
    {
        binprint(Base[i]);
        printf("\n");
    }
}

void ReadTargetMatrix()
{
    TheMatrix >> NumTargets;
    TheMatrix >> NumInputs;
    MaxDist = 0;

    if (NumInputs >= 8 * sizeof(long long int))
    {
        cout << "too many inputs" << endl;
        exit(0);
    }

    int bit;
    for (int i = 0; i < NumTargets; i++)
    {
        long long int PowerOfTwo = 1;
        Target[i] = 0;
        Dist[i] = -1;

        for (int j = 0; j < NumInputs; j++)
        {
            TheMatrix >> bit;
            if (bit)
            {
                Dist[i]++;
                Target[i] = Target[i] + PowerOfTwo;
            }
            PowerOfTwo = PowerOfTwo * 2;
        }
        if (Dist[i] > MaxDist)
            MaxDist = Dist[i];
        TheMatrix.get();
    }
}

bool is_target(long long int x)
{
    for (int i = 0; i < NumTargets; i++)
        if (x == Target[i])
            return true;
    return false;
}

bool is_base(long long int x)
{

    if (x == 0)
    {
        cout << "asking if 0 is in Base " << endl;
        exit(0);
    }

    for (int i = 0; i < BaseSize; i++)
        if (x == Base[i])
            return true;
    return false;
}

int NewDistance(int u)
{

    if (Target[u] == 0)
        return 0;
    else
    {
        if (is_base(Target[u]) || (NewBase == Target[u]))
            return 0;

        if (reachablE(Target[u] ^ NewBase, Dist[u] - 1, 0, pow(2, DepthLimit) - DepthNewBase))
            return (Dist[u] - 1);
        else
            return Dist[u];
    }
}

bool reachablE(long long int T, int K, int S, long long int L)
{
    if ((BaseSize - S) < K)
        return false;
    if (L < 1)
        return false;
    if (K == 0)
        return false;
    if (K == 1)
    {
        for (int i = S; i < BaseSize; i++)
            if ((T == Base[i]) && (pow(2, Depth[i]) <= L))
                return true;
        return false;
    }

    if (reachablE(T ^ Base[S], K - 1, S + 1, L - pow(2, Depth[S])))
        return true;

    if (reachablE(T, K, S + 1, L))
        return true;

    return false;
}
