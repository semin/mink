#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAXLINE 30000

double eucleadian_distance(int, double*, double*);
void chomp(char*);

int main(int argc, char *argv[]) {

    if (argc != 2) {
        printf( "Usage: %s descriptor_size < descriptor_file.csv\n", argv[0] );
        exit(1);
    }

    /* descriptor size */
    int ds = atoi(argv[1]);

    typedef struct {
        int no;
        double descriptor[ds];
    } Descriptor;

    Descriptor descriptors[MAXLINE];

    char line[MAXLINE];
    char* elem;
    char* delims = ", ";
    int dcount = 0;
    int p, q, r;
    double d;

    while(fgets(line, sizeof(line), stdin) != NULL) {
        Descriptor d;

        chomp(line);
        elem = strtok(line, delims);

        if (elem != NULL) {
            d.no = atoi(elem);
        } else {
            printf("Something wrong.\n");
            exit(1);
        }

        int j = 0;
        while (elem != NULL) {
            elem = strtok(NULL, delims);
            if ((elem != NULL) && (j < ds)) {
                d.descriptor[j] = atof(elem);
                j++;
            }
        }

        if (dcount < MAXLINE) {
            descriptors[dcount] = d;
            dcount++;
        } else {
            printf("You have more than %d lines!\n", MAXLINE);
        }
    }

    r = 0;
    for (p = 0; p < dcount - 1; p++) {
        for (q = p + 1; q < dcount; q++) {
            r++;
            d = eucleadian_distance(ds, descriptors[p].descriptor, descriptors[q].descriptor);
            printf("0,%d,%d,%f\n", descriptors[p].no, descriptors[q].no, d);
        }
    }

    return 0;
}

double eucleadian_distance(int ds, double* a, double* b) {
    int i = 0;
    double tmp_sum  = 0.0;
    double distance = 0.0;

    for(i = 0; i < ds; i++) {
        tmp_sum += pow(fabs(a[i] - b[i]), 2);
    }

    distance = sqrt(tmp_sum);

    return distance;
}

void chomp (char* s) {
    int end = strlen(s) - 1;
    if (end >= 0 && s[end] == '\n') {
        s[end] = '\0';
    }
}
