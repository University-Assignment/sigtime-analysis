from pdf2image import convert_from_path

pages = convert_from_path('C:/Users/ghks0/Desktop/Project/sigtime-server/uploads/60cecb7087b4342798c4da0a/results/SBS96/All_Solutions/SBS96_3_Signatures/Signatures/Signature_plotSBS_96_plots_S3.pdf', poppler_path="C:/Users/ghks0/Desktop/Project/sigtime-analysis/poppler-21.03.0/Library/bin")
i = 0
for page in pages:
    page.save('C:/Users/ghks0/Desktop/Project/sigtime-server/uploads/60cecb7087b4342798c4da0a/results/SBS96/All_Solutions/SBS96_3_Signatures/Signatures/Signature_plotSBS_96_plots_S3_%d.jpg' % i, 'JPEG')
    i += 1