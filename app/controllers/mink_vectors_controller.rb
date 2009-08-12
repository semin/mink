class MinkVectorsController < ApplicationController

  def index
    @mink_vectors = MinkVector.paginate(:page => params[:page] || 1, :per_page => 10)

    respond_to do |format|
      format.html
    end
  end

  def show
    @mink_vector = MinkVector.find(params[:id])

    respond_to do |format|
      format.html
    end
  end

  def search
    @query = params[:query]
    @mink_vector = MinkVector.search(@query,
                                     :match_mode => :extended,
                                     :page => params[:page],
                                     :per_page => 10)

    respond_to do |format|
      format.html
    end
  end

end
