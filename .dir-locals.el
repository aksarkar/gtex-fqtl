((org-mode . ((org-publish-project-alist .
                                         (("fqtl-org" . (:base-directory "/project2/mstephens/aksarkar/projects/gtex-fqtl/org"
                                                                         :publishing-directory "/project2/mstephens/aksarkar/projects/gtex-fqtl/docs"
                                                                         :publishing-function org-html-publish-to-html
                                                                         :exclude "setup.org"
                                                                         :htmlized-source t
                                                                         ))
                                          ("fqtl-fig" . (:base-directory "/project2/mstephens/aksarkar/projects/gtex-fqtl/org"
                                                                         :publishing-directory "/project2/mstephens/aksarkar/projects/gtex-fqtl/docs"
                                                                         :publishing-function org-publish-attachment
                                                                         :base-extension "png"
                                                                         :recursive t
                                                                         ))
                                          ("fqtl" . (:components ("fqtl-org" "fqtl-fig"))))))))
