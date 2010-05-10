class GitVectorSimilarity < ActiveRecord::Base

  belongs_to  :git_vector,
              :class_name   => "NormGitVector",
              :foreign_key  => "git_vector_id"

  belongs_to  :git_vector_target,
              :class_name   => "NormGitVector",
              :foreign_key  => "similar_git_vector_id"

end
