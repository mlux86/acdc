#include <unordered_map>

template<typename T>
Eigen::VectorXd Util::numericLabels(const std::vector<T> & labels)
{
    unsigned n = labels.size();
    Eigen::VectorXd v(n);
    std::unordered_map<T, unsigned> lblMap;
    unsigned cnt = 0;
    unsigned i = 0;
    for (const auto & lbl : labels)
    {
        if(lblMap.find(lbl) == lblMap.end())
        {
            lblMap[lbl] = cnt++;
        }
        v(i++) = lblMap[lbl];
    }   
    return v;	
}
