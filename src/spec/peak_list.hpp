#ifndef PROT_PEAK_LIST_HPP_
#define PROT_PEAK_LIST_HPP_

namespace prot {

template <class T>
bool comparePos(DeconvPeakPtr a, DeconvPeakPtr b);

template <class T>
T findHighestPeak(std::vector<T> ptr_list);

template <class T>
int findHighestPeakIdx(std::vector<T> ptr_list);

template <class T>
void sortOnPos(std::vector<T> &ptr_list);

template <class T>
double findMaxPos(std::vector<T> &ptr_list);

template <class T>
void sortOnIntensity(std::vector<T> &ptr_list);

template <class T>
int searchPos(std::vector<T> &ptr_list, double pos);

/**
 * Finds the highest peak for a specific position with error tolerance.
 */
template <class T>
int getHighPeakIdx(std::vector<T> &ptr_list, double pos, double tolerance);

/**
 * Finds the nearest peak for a specific position with error tolerance.
 */
template <class T>
int getNearPeakIdx(std::vector<T> &ptr_list, double pos, double tolerance);

/**
 * Removes all peaks in [center-interval, center+interval]
 */
template <class T>
std::vector<T> rmPeaks(std::vector<T> &ptr_list, double center, 
                       double interval);

template <class T>
std::vector<T> rmZeroPeaks(std::vector<T> &ptr_list);

/**
 * Removes a list of peaks.
 */
template <class T>
std::vector<T> rmPeaks(std::vector<T> &ptr_list, std::vector<bool> keep);

/**
 * Removes precursor mass. In top down data, one peak may be split to two
 * peaks. try to combine them together.
 */
template <class T>
std::vector<T> rmClosePeaks(std::vector<T> &ptr_list, double tolerance);
}

#endif
