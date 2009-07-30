class MinkVector < ActiveRecord::Base

  include Mink::Distance

  belongs_to  :scop_domain,
              :foreign_key  => :scop_id

  has_one :norm_mink_vector

end