function output = imbinarize2(image)
% adapted version of imbinarize to get rid of all-one image outputs
if graythresh(image)==0
    output = false(size(image));
else
    output = imbinarize(image);
end