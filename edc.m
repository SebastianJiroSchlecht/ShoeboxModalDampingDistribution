function e = edc(x,t)

e = db(flipud(cumsum(flipud(x.^2),1)) ./ sum(x(t:end,:).^2,1));