class ScopDomainsController < ApplicationController

  def index
    @scops = ScopDomain.representative.paginate(:page => params[:page] || 1,
                                                :per_page => 20)

    respond_to do |format|
      format.html
    end
  end

  def show
    @scop = ScopDomain.find_by_sunid(params[:id])

    respond_to do |format|
      format.html
    end
  end

  def search
    @query = params[:query]
    @scops = ScopDomain.representative.search(@query,
                                              :match_mode => :extended,
                                              :page => params[:page],
                                              :per_page => 10)

    respond_to do |format|
      format.html
    end
  end

end
