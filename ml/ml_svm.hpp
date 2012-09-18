#ifndef ML_SVM_HPP
#define ML_SYM_HPP

/*
Copyright (c) 2010 Fang-Cheng Yeh
All rights reserved.

Modified from LIBSVM Copyright (c) 2000-2010 Chih-Chung Chang and Chih-Jen Lin
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither name of copyright holders nor the names of its contributors
may be used to endorse or promote products derived from this software
without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include <memory>
#include <vector>

#define LIBSVM_VERSION 291

extern int libsvm_version;

struct svm_problem
{
	unsigned int l;
	unsigned int att_count;
	double *y;
	double **x;
};

enum { C_SVC, NU_SVC, ONE_CLASS, EPSILON_SVR, NU_SVR };	/* svm_type */
enum { LINEAR, POLY, RBF, SIGMOID, PRECOMPUTED }; /* kernel_type */

struct svm_parameter
{
	int svm_type;
	int kernel_type;
	int degree;	/* for poly */
	double gamma;	/* for poly/rbf/sigmoid */
	double coef0;	/* for poly/sigmoid */

	/* these are for training only */
	double cache_size; /* in MB */
	double eps;	/* stopping criteria */
	double C;	/* for C_SVC, EPSILON_SVR and NU_SVR */
	int nr_weight;		/* for C_SVC */
	int *weight_label;	/* for C_SVC */
	double* weight;		/* for C_SVC */
	double nu;	/* for NU_SVC, ONE_CLASS, and NU_SVR */
	double p;	/* for EPSILON_SVR */
	int shrinking;	/* use the shrinking heuristics */
	int probability; /* do probability estimates */
};

struct svm_model *svm_train(const struct svm_problem *prob, const struct svm_parameter *param);
void svm_cross_validation(const struct svm_problem *prob, const struct svm_parameter *param, int nr_fold, double *target);

int svm_get_svm_type(const struct svm_model *model);
int svm_get_nr_class(const struct svm_model *model);
void svm_get_labels(const struct svm_model *model, int *label);
double svm_get_svr_probability(const struct svm_model *model);

double svm_predict_values(const struct svm_model *model, const double *x, unsigned int att_count,double* dec_values);
double svm_predict(const struct svm_model *model, const double *x,unsigned int att_count);
double svm_predict_probability(const struct svm_model *model, const double *x, unsigned int att_count,double* prob_estimates);

void svm_destroy_model(struct svm_model *model);
void svm_destroy_param(struct svm_parameter *param);

const char *svm_check_parameter(const struct svm_problem *prob, const struct svm_parameter *param);
int svm_check_probability_model(const struct svm_model *model);

void svm_set_print_string_function(void (*print_func)(const char *));


template<typename attribute_type,typename classification_type>
class svm
{
	std::vector<double> y_buf;
	std::vector<std::vector<double> > data;
	std::vector<double*> x_buf;
	svm_parameter param;
	svm_problem prob;
	svm_model* model;
public:
	svm(void):model(0)
    {
        // default values
        param.svm_type = C_SVC;
        param.kernel_type = RBF;
        param.degree = 3;
        param.gamma = 0;	// 1/num_features
        param.coef0 = 0;
        param.nu = 0.5;
        param.cache_size = 100;
        param.C = 100;
        param.eps = 1e-3;
        param.p = 0.1;
        param.shrinking = 1;
        param.probability = 0;
        param.nr_weight = 0;
        param.weight_label = NULL;
        param.weight = NULL;
    }
	~svm(void)
	{
		svm_destroy_model(model);
	}
public:
    template<typename attributes_iterator_type,typename classifications_iterator_type>
    void learn(attributes_iterator_type attributes_from,
               attributes_iterator_type attributes_to,
               unsigned int attribute_dimension_,
               classifications_iterator_type classifications_from)
    {
        // put data in prob
		prob.att_count = attribute_dimension_;
        {
			prob.l = attributes_to-attributes_from;
            y_buf.resize(prob.l);
			x_buf.resize(prob.l);
			data.resize(prob.l);
            for (unsigned int index = 0;index < prob.l;++index)
            {
				data[index].resize(prob.att_count);
				std::copy(&attributes_from[index][0],&attributes_from[index][0]+prob.att_count,data[index].begin());
                x_buf[index] = &*(data[index].begin());
                y_buf[index] = classifications_from[index];
			}
			prob.y = &*y_buf.begin();
            prob.x = &*x_buf.begin();
        }
		param.gamma = 1.0/(double)prob.att_count;
		if(model)
			svm_destroy_model(model);
		model = svm_train(&prob,&param);
    }
    template<typename sample_iterator_type>
    classification_type predict(sample_iterator_type predict_attributes_) const
    {
		std::vector<double> x(predict_attributes_,predict_attributes_+prob.att_count);
        return svm_predict(model,&x[0],prob.att_count);
    }
};





#endif//ML_SYM_HPP
