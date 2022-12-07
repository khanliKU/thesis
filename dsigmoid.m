function r = dsigmoid(x)
    r = sigmoid(x).*(1-sigmoid(x));
end