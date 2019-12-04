function normonesmat = normonesmat(trap)

whts = whts(trap);
rnorms = normals(trap);
whts = whts(:);
whts2 = repmat(whts.',2,1);
whts2 = whts2(:).*rnorms(:);

normonesmat = bsxfun(@times,rnorms(:),whts2.');

end
