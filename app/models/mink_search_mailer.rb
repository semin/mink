class MinkSearchMailer < ActionMailer::Base

  def notification_email(mink_search)
    recipients  = mink_search.email
    from        = "MINK mailer <mink@cryst.bioc.cam.ac.uk>"
    subject     = "[MINK] Your job, #{mink_search.uuid} completed"
    sent_on     = Time.now
    body        :mink_search => mink_search
  end

end
