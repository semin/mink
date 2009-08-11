class MinkSearch < ActiveRecord::Base

  has_attached_file :pdb

  has_many :mink_search_hits

  def to_param
    self.uuid
  end
end
