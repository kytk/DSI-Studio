#ifndef ML_K_MEANS_HPP
#define ML_K_MEANS_HPP



template<typename attribute_type,typename classifications_iterator_type>
void k_means_clustering(const normalized_attributes<attribute_type>& attributes,
                        classifications_iterator_type classification,unsigned int k)
{
    unsigned int sample_size = attributes.size();
    unsigned int attribute_dimension = attributes.attribute_dimension();
    for (unsigned int index = 0;index < sample_size;++index)
        classification[index] = index%k;
    unsigned int change_cluster = 0;
    do
    {
        // E-step
        std::vector<std::vector<attribute_type> > means(k);
        for (unsigned int index = 0;index < k;++index)
        {
            std::vector<attribute_type> cur_mean(attribute_dimension);
            unsigned int count = 0;
            for (unsigned int j = 0;j < sample_size;++j)
                if (classification[j] == index)
                {
                    for (unsigned int k = 0;k < attribute_dimension;++k)
                        cur_mean[k] += attributes[j][k];
                    ++count;
                }
            if (count)
                for (unsigned int k = 0;k < attribute_dimension;++k)
                    cur_mean[k] /= (attribute_type)count;
            means[index].swap(cur_mean);
        }
        // M-step
        change_cluster = 0;
        for (unsigned int j = 0;j < sample_size;++j)
        {
            attribute_type min_dis = std::numeric_limits<attribute_type>::max();
            unsigned int min_cluster = 0;
            for (unsigned int index = 0;index < k;++index)
            {
                attribute_type dis2 = 0;
                for (unsigned int i = 0;i < attribute_dimension;++i)
                {
                    attribute_type d = attributes[j][i] - means[index][i];
                    dis2 += d*d;
                    if (dis2 > min_dis)
                        break;
                }
                if (dis2 < min_dis)
                {
                    min_dis = dis2;
                    min_cluster = index;
                }
            }
            if (classification[j] != min_cluster)
            {
                classification[j] = min_cluster;
                ++change_cluster;
            }
        }
    }
    while (change_cluster);
}


template<typename attribute_type,typename classification_type>
class k_means
{
protected:
    unsigned int k;
public:

public:
    k_means(unsigned int k_):k(k_) {}

    template<typename attributes_iterator_type,typename classifications_iterator_type>
    void operator()(attributes_iterator_type attributes_from,
                    attributes_iterator_type attributes_to,
                    unsigned int attribute_dimension,
                    classifications_iterator_type classifications_from)
    {
        normalized_attributes<attribute_type> attributes(attributes_from,attributes_to,attribute_dimension);
        k_means_clustering(attributes,classifications_from,k);
    }

    template<typename attributes_iterator_type,typename classifications_iterator_type>
    void operator()(attributes_iterator_type attributes_from,
                    attributes_iterator_type attributes_to,
                    classifications_iterator_type classifications_from)
    {
        normalized_attributes<attribute_type> attributes(attributes_from,attributes_to);
        k_means_clustering(attributes,classifications_from,k);
    }
};

#endif//ML_K_MEANS_HPP
