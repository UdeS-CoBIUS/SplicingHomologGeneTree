#ifndef MEMORY_H_INCLUDED
#define MEMORY_H_INCLUDED


short int ** allocateTwoDim_short(unsigned int row, unsigned int col);
void destroyTwoDim_short(short ** ptr, unsigned row, unsigned int col);
float ** allocateTwoDim_float(unsigned row, unsigned col);
void destroyTwoDim_float(float ** ptr, unsigned row, unsigned col);


#endif // memory_H_INCLUDED
