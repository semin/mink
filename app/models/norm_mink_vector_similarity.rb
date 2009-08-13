class NormMinkVectorSimilarity < ActiveRecord::Base

  belongs_to  :norm_mink_vector,
              :class_name   => "NormMinkVector",
              :foreign_key  => "norm_mink_vector_id"

  belongs_to  :norm_mink_vector_target,
              :class_name   => "NormMinkVector",
              :foreign_key  => "similar_norm_mink_vector_id"

end
