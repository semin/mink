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
    @mink_vector = MinkVector.find(params[:id])
    @scop_domain = @mink_vector.scop_domain

    respond_to do |format|
      format.html
    end
  end

end
