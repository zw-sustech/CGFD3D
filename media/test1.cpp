#include<iostream>
#include<cmath>
//#include<cfloat>
using namespace std;

struct tt{
    int *a;
    int *b;
};

void p(int *a, int n) {
    for (int i = 0; i < n; i++) {
        cout << *(a+i) << endl;
//        scanf("%d",a+i);
    }
}

void init(int slice, int ni, tt *t) {
    t->a = new int[slice * ni];
    t->b = new int[slice * ni];
}

void write(int *u, int slice, int val) {
    for (int i = 0; i < slice; i++)
        scanf("%d",u+i);
}

int main()
{
    tt t;
    int NI = 3, slice = 6;

    init(slice, NI, &t);

    for (int i = 0; i < NI; i++) {
        write(t.a+slice*i, slice, i);
    }

    int ni = 2;
    p(t.a+ni*slice, slice);
    p(t.b+ni*slice, slice);
    return 0;
}
