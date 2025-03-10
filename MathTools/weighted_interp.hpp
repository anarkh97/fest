using namespace std;

namespace MathTools {

// helper functions
double vec_distance(int n, double v1[], double v2[]);

//interpolation functions
void weighted_interp(int dim, int nd, double xd[], double fd[],
  int nq, double xq[], double fq[]);

}
