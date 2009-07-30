class NormMinkVector < ActiveRecord::Base

  include Mink::Distance

  belongs_to  :mink_vector

end
