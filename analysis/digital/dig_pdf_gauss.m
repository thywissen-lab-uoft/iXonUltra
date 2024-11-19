function digdata = dig_pdf_gauss(digdata,opts)

Xc_pdf = zeros(size(digdata.Zdig,3),2);
Xs_pdf = zeros(size(digdata.Zdig,3),2);
Yc_pdf = zeros(size(digdata.Zdig,3),2);
Ys_pdf = zeros(size(digdata.Zdig,3),2);
for kk=1:size(digdata.Zdig,3)
    Ratom = digdata.Ratom{kk}; % Position of very atom
    X = Ratom(1,:);
    Y = Ratom(2,:);
          
    [pdf_x,pdf_x_cint] = mle(X,'distribution','normal');
    [pdf_y,pdf_y_cint] = mle(Y,'distribution','normal');

    Xc_pdf(kk,1)=   pdf_x(1);
    Xc_pdf(kk,2)=  0.5*(pdf_x_cint(2,1)-pdf_x_cint(1,1));
    Xs_pdf(kk,1)=   pdf_x(2);
    Xs_pdf(kk,2)=  0.5*(pdf_x_cint(2,2)-pdf_x_cint(1,2));

    Yc_pdf(kk,1)=   pdf_y(1);
    Yc_pdf(kk,2)=  0.5*(pdf_y_cint(2,1)-pdf_y_cint(1,1));
    Ys_pdf(kk,1)=   pdf_y(2);
    Ys_pdf(kk,2)=  0.5*(pdf_y_cint(2,2)-pdf_y_cint(1,2));      
end
digdata.Xc_pdf_px = Xc_pdf;
digdata.Yc_pdf_px = Yc_pdf;
digdata.Xs_pdf_px = Xs_pdf;
digdata.Ys_pdf_px = Ys_pdf;
end