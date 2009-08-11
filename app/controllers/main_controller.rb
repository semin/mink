class MainController < ApplicationController

  def home
    @recent_news = News.order('date desc').limit(5)

    respond_to do |format|
      format.html
    end
  end

  def references
    respond_to do |format|
      format.html
    end
  end

  def contacts
    respond_to do |format|
      format.html
    end
  end

  def news
    @news = News.order('date desc')
  end

  def search
    redirect_to :controller => "scop_domains",
                :action => "search",
                :query => params[:query]
  end
end
