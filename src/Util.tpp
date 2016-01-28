#include <unordered_map>

template<typename T>
std::vector<unsigned> Util::numericLabels(const std::vector<T> & labels)
{
    unsigned n = labels.size();
    std::vector<unsigned> v(n, 0);
    std::unordered_map<T, unsigned> lblMap;
    unsigned cnt = 0;
    unsigned i = 0;
    for (const auto & lbl : labels)
    {
        if(lblMap.find(lbl) == lblMap.end())
        {
            lblMap[lbl] = cnt++;
        }
        v[i++] = lblMap[lbl];
    }   
    return v;	
}

template<typename T>
std::vector<T> Util::mutualLabels(const std::vector<T> & labels1, const std::vector<T> & labels2)
{
    auto l1 = labels1;
    auto l2 = labels2;

    std::sort(l1.begin(), l1.end());
    std::sort(l2.begin(), l2.end());

    std::vector<T> v(std::max(labels1.size(), labels2.size()));
    auto it = std::set_intersection(l1.begin(), l1.end(), l2.begin(), l2.end(), v.begin());
    v.resize(it - v.begin());

    return v;
}