module Matrix where

type Matrix = (Int, [Float])
zeroes     = (iterate (+0) 0)

-- Making matrices

-- Creates a row matrix from a list of floats
-- Fills the rest with 0
createMat :: [Float] -> Matrix
createMat xs
    = (size, full)
    where
        size       = length xs
        full       = xs ++ take (size * (size - 1)) zeroes

-- Creates a column matrix from a list of floats
createColMat :: [Float] -> Matrix
createColMat = transposeM . createMat

--Creates a list from a column matrix
createListfromColMatrix :: Matrix -> [Float]
createListfromColMatrix m@(s, xs)
    = take s ((snd . transposeM) m)

colList :: Matrix -> [Float]
colList = createListfromColMatrix

-- Useful utility functions for matrices 

-- Applies a function to every element in the matrix
mapMat :: (Float -> Float) -> Matrix -> Matrix
mapMat f (s, xs) = (s, map f xs)

-- Makes the matrix smaller if the outermost layers are all zero
decrease :: Matrix -> Matrix
decrease m@(s, xs)
    | s == 1       = m
    | lrow && lcol = decrease (reduce m (s - 1) (s - 1))
    | otherwise    = m
        where
            (_, ts) = transposeM m
            lrow = eqZero ((split s xs) !! (s - 1))
            lcol = eqZero ((split s ts) !! (s - 1))
            eqZero  = and . (map (==0))

-- Takes a list and splits into lists of size n
split :: Int -> [Float] -> [[Float]]
split n xs 
    | null xs   = []
    | otherwise =  y : split n y'
        where
            (y, y') = splitAt n xs

-- Takes a matrix and returns list of rows of the matrix
splitMat :: Matrix -> [[Float]]
splitMat (n, xs) = split n xs

-- Transposes a list
transpose :: [[Float]] -> [[Float]]
transpose xs
    = [map (!!i) xs | i <- [0.. length xs - 1]]

-- Transposes a matrix (flip elements across diagonal)
transposeM :: Matrix -> Matrix
transposeM m@(n, _)
    = (n, concat (transpose (splitMat m)))

-- Deletes all the elements on given row and column
reduce :: Matrix -> Int -> Int -> Matrix
reduce (m, xs) i j 
    = (m - 1, reduced)
        where
            reduced = [xs !! n | n <- [0..m * m - 1], let (q, r) = quotRem n m, q /= i, r /= j]

-- Basic Matrix operations

--Adds 2 matrices of same size together
addMatrix :: Matrix -> Matrix -> Matrix
addMatrix (m1, m1M) (m2, m2M)
    | m1 /= m2  = error "Matrix dimensions differ!"
    | otherwise = (m1, zipWith (+) m1M m2M)

-- Scalar multiplication on matrix with given scalar
scalarMult :: Float -> Matrix -> Matrix
scalarMult n
    = mapMat (*n)

-- Multiplies two matrices if they have same dimensions
multMatrix :: Matrix -> Matrix -> Matrix
multMatrix a@(m1, _) b@(m2, _) 
    | m1 /= m2  = error "Matrix dimensions differ!"
    | otherwise = multMatrix' a b
    where
        multMatrix' :: Matrix -> Matrix -> Matrix
        multMatrix' (m, m1M) (_, m2M)
            = (m, res)
            where
                mat1 = split m m1M
                mat2 = transpose (split m m2M)
                res = [sum (zipWith (*) is js) | is <- mat1, js <- mat2]

-- Infix notation for above operations

infixr 8 ^^^
(^^^) = powerMatrix

infixl 7 ***
(***) = multMatrix

infixl 7 \**
(\**) = scalarMult

infixl 7 +++
(+++) = addMatrix

-- Returns the identity matrix of dimension n
identity :: Int -> Matrix
identity n 
    = (n, [if (q == r) then 1 else 0 | i <- [0..n*n - 1], let (q, r) = quotRem i n])

-- Returns the matrix calculated by multiplying the matrix m k times
powerMatrix :: Matrix -> Int -> Matrix
powerMatrix m@(n, _) k 
    = foldl multMatrix (identity n) [m | x <- [1..k]]

-- Returns the determinant of a matrix
determinant :: Matrix -> Float
determinant (1, x) = head x
determinant m@(n, xs)
    = sum (zipWith (*) coeffs reducedDets)
        where
            signs = [(-1) ^ (x + 1) | x <- [1..]]
            coeffs = zipWith (*) (head (splitMat m)) signs
            reducedM = [reduce m 0 j | j <- [0..n-1]]
            reducedDets = map determinant reducedM

-- Inverts a matrix if determinant is non-zero
invert :: Matrix -> Matrix
invert m@(n, _)
    | determinant m == 0 = error "Matrix has no inverse!"
    | otherwise          = coeff \** mTrans
        where
            minors  = [determinant (reduce m i j) | i <- [0..n-1], j <- [0..n-1]]
            mCofact = zipWith (*) minors [(-1) ^ (x + y) | x <- [1..], y <- [1..n]]
            mTrans  = transposeM (n, mCofact)
            coeff   = (1 / determinant m)

--Returns the transitive closure of a matrix
transClose :: Matrix -> Matrix
transClose m@(n, _)
    = foldl1 addMatrix as
        where
            as = take n (iterate (*** m) m)

-- Returns a matrix of functions
funcMatrix ::  Int -> [Float -> Float] -> Float -> Matrix
funcMatrix n fs x
    | length fs /= n*n = error "Size of list does not match"
    | otherwise        = (n, res)
    where
        res = map (\f -> f x) fs

-- Rotation matrix of size n of angle θ
rotMatr2 :: Float -> Matrix
rotMatr2 θ 
    = funcMatrix 2 [cos, negate . sin, sin, cos] θ

-- Rotates a point in ℝ² by angle θ
rot2 :: (Float, Float) -> Float -> (Float, Float)
rot2 (x, y) θ
    = (x', y')
    where
        point = createColMat [x, y]
        rot = rotMatr2 θ
        [x', y'] = createListfromColMatrix (rot *** point)

-- Increases matrix of size n to n+1
inc1 :: Matrix -> Matrix
inc1 m@(n, xs)
    = (n + 1, concat xs')
    where
        addZeroes = map (\x -> x ++ [0]) (splitMat m)
        xs' = addZeroes ++ [take (n+1) zeroes]

-- Increases matrix of size n to size m
inc :: Matrix -> Int -> Matrix
inc mat@(n, _) m
    | m == n    = mat
    | n > m     = error "Matrix larger than given size"
    | otherwise = inc (inc1 mat) m


-- Strassen's Algorithm for Matrix multiplication
strassMult2 :: Matrix -> Matrix -> Matrix
strassMult2 a@(n1, as) b@(n2, bs)
    | n1 /= 2 || n2 /= 2 = error "Size not valid"
    | otherwise          = (2, [x1 + x4 -x5 + x7, x3 + x5, x2 + x4, x1 + x3 - x2 + x6])
    where
        x1 = (as !! 0 + as !! 3) * (bs !! 0 + bs !! 3)
        x2 = (as !! 2 + as !! 3) * (bs !! 0)
        x3 = (as !! 0) * (bs !! 1 - bs !! 3)
        x4 = (as !! 3) * (bs !! 2 - bs !! 0)
        x5 = (as !! 0 + as !! 1) * (bs !! 3)
        x6 = (as !! 2 - as !! 0) * (bs !! 0 + bs !! 1)
        x7 = (as !! 1 - as !! 3) * (bs !! 2 + bs !! 3)

-- strassMult :: Matrix -> Matrix -> Matrix
-- strassMult m1@(n1, _) m2@(n2, _)
--     | n1 /= n2 = error "Matrix not of same size"
--     | otherwise = decrease (strassMult' (inc m1 n) (inc m2 n))
--         where
--             n = head [ j | i <- [1..], let j = 2^i, j >= n1]
--             strassMult' :: Matrix -> Matrix -> Matrix
--             strassMult a@(2, _) b 
--                 = strassMult2 a b
--             strassMult a@(n, as) bs(_, bs)
--                 = 
--                     where
--                         [a1,a2,a3,4] = split (n/2) as
--                         bs' = split (n/2) bs