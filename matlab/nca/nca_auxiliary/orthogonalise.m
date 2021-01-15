function v_ort = orthogonalise(v,W)

v_ort = v-(W'*pinv(W')*v')';
