class MinkSearchesController < ApplicationController

  before_filter :authenticate, :only => [ :index, :destroy ]

  def index
    @mink_searches = MinkSearch.all

    respond_to do |format|
      format.html # index.rhtml
    end
  end

  def show
    if @mink_search = MinkSearch.find_by_uuid(params[:id])
      @sorted_mink_search_hits = if !@mink_search.finished_at.nil?
                                  @mink_search.mink_search_hits.sort_by(&:distance).slice(0, 50)
                                end

      respond_to do |format|
        format.html
      end
    else
      redirect_to "/"
    end
  end

  def new
    @mink_search = MinkSearch.new
  end

  def create
    @mink_search        = MinkSearch.new(params[:mink_search])
    @mink_search.uuid   = UUID.new.generate(:compact)
    @mink_search.status = "Standing by"

    respond_to do |format|
      if @mink_search.save
        flash[:notice] = 'Your job was succesfully submitted.'
        Delayed::Job.enqueue MinkSearchJob.new(@mink_search.id)
        format.html { redirect_to mink_search_url(@mink_search.uuid) }
      else
        format.html { render :action => 'new' }
      end
    end
  end

  def destroy
    @mink_search = MinkSearch.find_by_uuid(params[:id])
    @mink_search.destroy

    respond_to do |format|
      format.html { redirect_to mink_searches_url }
      format.xml { render :nothing => true }
    end
  end

  protected

  def authenticate
    authenticate_or_request_with_http_basic do |username, password|
      Digest::SHA1.hexdigest(username) == configatron.username &&
        Digest::SHA1.hexdigest(password) == configatron.password
    end
  end
end
