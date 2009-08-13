class MinkVectorsController < ApplicationController

  def index
    @query = params[:query]

    if @query && !@query.empty?
      @mink_vectors = MinkVector.search(@query,
                                        :match_mode => :extended,
                                        :page => params[:page] || 1,
                                        :per_page => 10)
    else
      @mink_vectors = MinkVector.paginate(:page => params[:page] || 1,
                                          :per_page => 10)
    end

    respond_to do |format|
      format.html
    end
  end

  def show
    @mink_vector      = MinkVector.find(params[:id])
    @norm_mink_vector = @mink_vector.norm_mink_vector
    @scop_domain      = @mink_vector.scop_domain
    @sorted_similar_norm_mink_vectors = @mink_vector.
                                        norm_mink_vector.
                                        sorted_similar_norm_mink_vectors.slice(0, 50)

    respond_to do |format|
      format.html
    end
  end

end
