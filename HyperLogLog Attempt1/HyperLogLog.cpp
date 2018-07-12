
#include "HyperLoglog.h"
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include "MurMurHash3.h"

#include <math.h>
#include <omp.h>


HyperLogLog::HyperLogLog(unsigned int _bitsTaken) {
	bitsTaken = _bitsTaken;
	cardinality = 0;
}




void HyperLogLog::Add(char param[], int dna_sequence_size) {
	int numThreads = 0;
	numBuckets = pow(2, bitsTaken);

#pragma omp parallel
	{
		int threadNo = omp_get_thread_num();
#pragma omp single 
		{
			numThreads = omp_get_num_threads();
			buckets = new int*[numThreads];
		}

		buckets[threadNo] = new int[numBuckets];
		for (int i = 0; i < numBuckets; i++) {
			buckets[threadNo][i] = 0;
		}

#pragma omp for schedule(static)
		for (int i = 0; i < dna_sequence_size - windowSize + 1; ++i)
		{
			char* p;
			p = &param[i];
			//unsigned int encoded = Encode(p);
			//printf("%d\n", new_group);//
			unsigned int hash;
			MurmurHash3_x86_32(p, windowSize, hashSeed, (void*)&hash);
			unsigned int index = hash >> (32 - bitsTaken);
			//unsigned int buckets;
			int z = 0;
			int counter = 0;
			int size = (sizeof(unsigned int) * 8);
			unsigned int temp2 = (1 << 31);
			while (z < size) {
				unsigned int temp = (hash << (bitsTaken));

				if ((temp&temp2) == 0) {
					counter++;
				}
				else if ((temp&temp2) != 0) {

					z = size;

				}
				temp2 = (temp2 >> 1);
				z++;
			}

			if (buckets[threadNo][index] < counter + 1) {
				buckets[threadNo][index] = counter + 1;
			}
		}
	}
	int len = sizeof(buckets);
	for (int i = 1; i < len; i++)
	{
		for (int j = 0; j < numBuckets; j++)
		{
			if (buckets[0][j] < buckets[i][j]) {
				buckets[0][j] = buckets[i][j];
			}
		}
	}
}

void HyperLogLog::Print() {
	int i = 0;
	while (i < numBuckets) {
		cout << buckets[0][i] << " ";
		i++;
	}

	cout << endl << cardinality;
}

void HyperLogLog::EstimationEQ() {
	{
		double sum = 0;
		double estimation;
		for (int i = 0; i < numBuckets; i++)
		{
			sum = sum + pow(2, -buckets[0][i]);
		}

		double alpha;
		switch (numBuckets) {
		case 16:
			alpha = 0.673;
			break;
		case 32:
			alpha = 0.697;
			break;
		case 64:
			alpha = 0.709;
			break;
		default:
			alpha = 0.7213 / (1.0 + 1.079 / numBuckets);
			break;
		}
		estimation = alpha*numBuckets*numBuckets*(1 / sum);
		cardinality += ceil(estimation);

	}
}



unsigned int HyperLogLog::Encode(char * p) {
	unsigned int encoded;

	for (uint32_t i = 0; i < windowSize; ++i)
	{

		//printf("%c\n", (*p));
		//printf("%d\n", (*p) >> 1 & 3);
		encoded ^= ((*p++) >> 1 & 3) << (windowSize * 2 - 2 * (i + 1));

		/*
		--> (*p) >> 1 & 3) does the following:
		A 65 01000|00|1  0
		C 67 01000|01|1  1
		G 71 01000|11|1  3
		T 84 01010|10|0  2

		--> (K * 2 - 2 * (i + 1)) does the following:
		i.e. streamed input is "AGCG" (K=4), you'd want to shift 00(A) value 6 bits left, 11(G) value 4 bits left etc.
		*/
	}
	return encoded;

}