#include "spline.h"

#include <numeric> // for std::accumulate

namespace towr{

Spline::Spline(const VecTimes& poly_durations, int n_dim)
{
    uint n_polys = poly_durations.size();

    // 把n个三次样条曲线赋值给调用者
    cubic_polys_.assign(n_polys, CubicHermitePolynomial(n_dim));
    
    for(uint i=0; i<cubic_polys_.size(); ++i) {
        cubic_polys_.at(i).SetDuration(poly_durations.at(i));
    }

    UpdatePolynomialCoeff();
}

int Spline::GetSegmentID(double t_global, const VecTimes& durations)
{
    double eps = 1e-10; // double precision
    assert(t_global >= 0.0);

    double t = 0;
    int i=0;
    
    for(double d: durations) {
        
        t += d;

        if(t >= t_global-eps) // at junctions, returns previous spline (=)
            return i;

        i++;
   }

   assert(false); // this should never be reached
}

std::pair<int, double> Spline::GetLocalTime(double t_global, const VecTimes& durations) const
{
    
    int id = GetSegmentID(t_global, durations);

    double t_local = t_global;
    for(int i=0; i<id; i++)
        
        t_local -= durations.at(i);

    return std::make_pair(id, t_local);
}

const State Spline::GetPoint(double t_global) const
{
    int id; 
    double t_local;
    
    std::tie(id, t_local) = GetLocalTime(t_global, GetPolyDurations());

    return GetPoint(id, t_local);
}

const State Spline::GetPoint(int poly_id, double t_local) const
{
    return cubic_polys_.at(poly_id).GetPoint(t_local);
}

void Spline::UpdatePolynomialCoeff()
{
    for(auto& p : cubic_polys_)
        p.UpdateCoeff();
}

int Spline::GetPolynomialCount() const
{
    return cubic_polys_.size();
}

Spline::VecTimes Spline::GetPolyDurations() const
{
    VecTimes poly_durations;
    for (const auto& p : cubic_polys_)
        poly_durations.push_back(p.GetDuration());

    return poly_durations;
}

double Spline::GetTotalTime() const
{
  auto v = GetPolyDurations();
  
  return std::accumulate(v.begin(), v.end(), 0.0);
}

} /* namespace towr */

