class MinkVectorSimilarity < ActiveRecord::Base

  belongs_to  :mink_vector,
              :class_name   => "NormMinkVector",
              :foreign_key  => "mink_vector_id"

  belongs_to  :mink_vector_target,
              :class_name   => "NormMinkVector",
              :foreign_key  => "similar_mink_vector_id"

end
