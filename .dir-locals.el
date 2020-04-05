((org-mode . ((org-html-htmlize-output-type . "css")
              (org-publish-project-alist .
                                         (("fqtl-org" . (:base-directory "/broad/compbio/aksarkar/projects/gtex-fqtl/org"
                                                                         :publishing-directory "/broad/compbio/aksarkar/projects/gtex-fqtl/docs"
                                                                         :publishing-function org-html-publish-to-html
                                                                         :exclude "setup.org"
                                                                         :htmlized-source t
                                                                         ))
                                          ("fqtl-fig" . (:base-directory "/broad/compbio/aksarkar/projects/gtex-fqtl/org"
                                                                         :publishing-directory "/broad/compbio/aksarkar/projects/gtex-fqtl/docs"
                                                                         :publishing-function org-publish-attachment
                                                                         :base-extension "png"
                                                                         :recursive t
                                                                         ))
                                          ("fqtl" . (:components ("fqtl-org" "fqtl-fig")))))
              ))
 (python-mode . ((python-indent-offset 2))))
