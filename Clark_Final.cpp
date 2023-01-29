// Final Project
// Autumn Clark
// Sample from The Watt Spectrum

#include <cmath>
#include <iostream>
#include <functional>
#include <random> // For random numbers
#include <fstream> // Allows reading and writing to files
#include <boost/math/tools/roots.hpp> // has bisect, newton
//#include <boost/math/differentiation/finite_difference.hpp>
namespace tools = boost::math::tools;

// Code for e found at https://stackoverflow.com/questions/18773343/how-to-calculate-euler-constant-or-euler-powered-in-c/18773495
const double EulerConstant = std::exp(1.0);

// Method that returns the approximate value of a specified input E
//      This approximate distribution is used to approximate the Watt Distribution
double f(double E){
  // 0.4865sinh(sqrt(2E))e^(-E)
  //return E*E;
  return 0.4865*std::sinh(sqrt(2*E))*std::pow(EulerConstant, (-1*E));
}

// Method that returns the derivative value at E based on finite difference approximation
double FuncDeriv(double E){
  // I attempted to use the boost function, but the header file required doesn't work
  //      so I wrote my own numerical approximation
  //return boost::math::differentiation::finite_difference_derivative(f, E);

  double h = 0.00000000000001; // h should be very small
  return ((f(E+h)-f(E))/h);
}

// Termination function for bisect
bool TerminationFunc(double min, double max){
  return abs(max-min) <= 0.00001; // 1e-5 precision
}

// Method that performs numerical intergration by Simpson's Rule
//     a and b are the start and end of the interval
//     n is the number of bins to be used -- directly tied to accuracy -- must be even for Simpson's
//     function Function() is passed as an argument
double SimpsonInteg(double a, double b, int n, std::function<double(double x)> f){
  if(n%2 != 0){ // Check if n is odd
    std::cout << "Bins must be even to use Simpson's Rule!" << std::endl;
    return -1; // If odd, return -1
  }
  double delta = (b-a)/n; // Width of each trapeziod
  double total = f(a); // First call to function
  for(int i=1; i<n; i++){
    if(i%2 == 0){ // Even calls are multiplied by 2
      total += 2 * f(a + delta*(i));
    } else{ // Odd calls are multiplied by 4
      total += 4 * f(a + delta*(i));
    }
  }
  total += f(b); // Last call to the function
  total = (total * (delta/3)); // Multiply running total by delta / 3, then reassign
  return total;
}

// Method that writes bin info to a file
void WriteBins(std::vector<int> const& bins){
  std::ofstream file; // Assigns the filestream as file
  file.open("BinData.txt"); // Open the appropiate file for number count
  file << "Bin\tCount\n"; // Add the column labels
  for(int i = 0; i < bins.size(); i++){
    // Loop through all the bins and write them to the file
    file << i << "\t" << bins.at(i) << "\n";
  }
  file.close(); // Close the file
}

int main(){
  /* // Test code
  double E = 5; double Y = f(E);
  std::cout << "The Watt Distribution of " << E << " is " << Y << std::endl;
  std::cout << "Euler Constant: " << EulerConstant << std::endl;
  */
  // Part 1 - b //
  std::cout << "Simpson's Rule Integration over [0,10]\n\twith 100 bins: " << SimpsonInteg(0,10,100,f) << std::endl;

  // Part 2 //
  std::pair<double, double> range = tools::bisect(FuncDeriv, 0.0, 8.0, TerminationFunc);
  std::cout << "The range by bisect method: " << range.first << " , " << range.second << std::endl;
  double root = range.first + (range.second - range.first)/2;
  double max = f(root);
  std::cout << "The root by bisect method: " << root << std::endl;
  std::cout << "Maximum of Watt Distribution: " << max << std::endl;

  // Parts 3, 4, and 5 //
  std::vector<double> sample; // Vector to hold the y values of the sample points
  std::random_device rd; // Random numbers based on time for random seed
  std::default_random_engine generator(rd()); // Define the generator for random numbers
  std::uniform_real_distribution<double> distx(0,8); // Specify the distribution of random numbers for x values
  std::uniform_real_distribution<double> disty(0,max); // Specify the distribution of random numbers for y values

  while(sample.size() != 100000){ // Loop until sample has 100,000 values
    //create point
    std::pair<double,double> point;
    // randomly choose x and y (point first and second)
    point.first = distx(generator); // Returns a random number between 0 and 8
    point.second = disty(generator); // Returns a random number between 0 and the max of f
    if(point.second <= f(point.first)){ // Check if the point is under the curve
      sample.push_back(point.first); // If so, add to sample
    }
  }
  std::cout << "Size of sample: " << sample.size() << std::endl;

  // Part 6 //
  std::vector<int> bins; // Vector to hold the count of each bin
  for(int i=0; i<100; i++){ // Add 100 bins
    bins.push_back(0); // All counts start at 0
  }
  std::sort(begin(sample), end(sample)); // Sort samples so that they can be placed in order
  /* // Test code to look at sample values
  for(int i=95000; i<100000; i++){
    std::cout << "Sample " << i << " : " << sample.at(i) << std::endl; // Test code
  }
  */
  int index = 0; // Index for the current sample being binned
  int bin_multiplier = 1; // Multipler used to determine which bin we are looking at

  while(index < sample.size()){ // Loop until all samples have been binned
    double bin_end = ((bin_multiplier)*(0.08)); // Bin width is 0.08 ((8-0)/100)
    if(sample.at(index) <= bin_end){ // If sample is less than or equal to the value of bin_end
      bins.at(bin_multiplier-1) = (bins.at(bin_multiplier-1))+1; // Increment count of that bin
      index++; // Look at next sample
    } else{ // Doesn't go in this bin
      bin_multiplier++; // Look at next bin
    }
  }
  // Loop through all bins and display their counts
  for(int i=0; i<bins.size(); i++){
    std::cout << "Bin " << i << " : " << bins.at(i) << std::endl; // Test code
  }

  // Part 7 //
  WriteBins(bins);

  return 0;
}
